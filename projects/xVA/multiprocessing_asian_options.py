import numpy as np
import math
import matplotlib.pyplot as plt
import time 
from scipy.stats import norm
from multiprocessing import get_context
from concurrent.futures import ProcessPoolExecutor, as_completed

def gbm_path_sim_vectorized(
    s_start, r, sigma, T_sim, n_sims, m_steps, Z_matrix=None, rng=None
):
    """
    Vectorized GBM.
    """
    dt = T_sim / m_steps

    # RNG
    if Z_matrix is None:
        if rng is None:
            Z = np.random.normal(0.0, 1.0, (n_sims, m_steps))
        else:
            Z = rng.normal(0.0, 1.0, size=(n_sims, m_steps))
    else:
        if Z_matrix.shape != (n_sims, m_steps):
            raise ValueError(
                f"Z_matrix shape {Z_matrix.shape} does not match (n_sims, m_steps) {(n_sims, m_steps)}"
            )
        Z = Z_matrix

    log_returns = (r - 0.5 * sigma**2) * dt + sigma * math.sqrt(dt) * Z
    cumulative_log_returns = np.cumsum(log_returns, axis=1)

    simulated_paths = np.zeros((n_sims, m_steps + 1))
    simulated_paths[:, 0] = s_start
    simulated_paths[:, 1:] = s_start * np.exp(cumulative_log_returns)
    return simulated_paths


def price_asian_option_mc(
    s_current, k, r, sigma, T_maturity, n_sims, m_total_steps, omega,
    t_current=0.0, historical_prices=None, Z_matrix=None, rng=None,
    sum_known_override=None, m_known_steps_override=None
):
    """
    Monte Carlo Asian option pricer for multiprocessing.
    New: rng for reproducibility, and optional (sum_known, m_known_steps) overrides to avoid repeated slicing/summing.
    """

    if omega not in (1, -1):
        raise ValueError("omega must be 1 (call) or -1 (put)")

    # Determine state (t=0 or t>0)
    if (sum_known_override is not None) and (m_known_steps_override is not None):
        sum_known = float(sum_known_override)
        m_known_steps = int(m_known_steps_override)
        T_remaining = T_maturity - t_current
    else:
        if historical_prices is None or len(historical_prices) == 0:
            sum_known = 0.0
            m_known_steps = 0
            T_remaining = T_maturity - t_current
        else:
            hp = np.asarray(historical_prices)
            sum_known = float(np.sum(hp))
            m_known_steps = int(len(hp))
            T_remaining = T_maturity - t_current

    m_remaining_steps = m_total_steps - m_known_steps

    # Simulate futures if needed
    if m_remaining_steps > 0:
        all_future_paths = gbm_path_sim_vectorized(
            s_current, r, sigma, T_remaining, n_sims, m_remaining_steps,
            Z_matrix=Z_matrix, rng=rng
        )
        simulated_future_prices = all_future_paths[:, 1:]
        sum_future = np.sum(simulated_future_prices, axis=1)
        total_sum = sum_known + sum_future
        s_ave = total_sum / m_total_steps
    else:
        # Average already known
        s_ave = sum_known / m_total_steps

    payoffs = np.maximum(omega * (s_ave - k), 0)
    price = np.mean(payoffs) * np.exp(-r * T_remaining)
    return price


# Worker to price a chunk of paths at a fixed time j
def _price_chunk_for_time_light(
    j, t_current,
    s_current_chunk,          # shape: (chunk_len,)
    sum_known_chunk,          # shape: (chunk_len,)
    m_known_steps,
    k, r, sigma, T_maturity, m_steps, omega, n_sims_inner,
    i_start, base_seed
):
    chunk_len = s_current_chunk.shape[0]
    out = np.empty(chunk_len, dtype=float)

    for local_idx in range(chunk_len):
        i = i_start + local_idx
        rng = np.random.default_rng(np.random.SeedSequence([base_seed, j, i]))

        out[local_idx] = price_asian_option_mc(
            s_current=s_current_chunk[local_idx],
            k=k, r=r, sigma=sigma,
            T_maturity=T_maturity,
            n_sims=n_sims_inner,
            m_total_steps=m_steps,
            omega=omega,
            t_current=t_current,
            historical_prices=None,
            rng=rng,
            sum_known_override=sum_known_chunk[local_idx],
            m_known_steps_override=m_known_steps
        )
    return i_start, out


# Parallel driver 
def calculate_mtm_matrix_asian_parallel(
    all_paths: np.ndarray,
    k: float,
    r: float,
    sigma: float,
    T_maturity: float,
    omega: int,
    n_sims_inner: int,
    max_workers: int | None = None,
    chunk_size: int = 1000,
    base_seed: int = 42,
    verbose: bool = True,
):
    n_sims_outer, m_steps_plus_1 = all_paths.shape
    m_steps = m_steps_plus_1 - 1
    time_vector = np.linspace(0, T_maturity, m_steps + 1)
    mtm_matrix = np.zeros_like(all_paths)

    # Precompute prefix sums of prices S[:,1:] along time once
    # cumsum_S[i, j-1] = sum of S_i at times 1..j  (so j -> index j-1)
    cumsum_S = np.cumsum(all_paths[:, 1:], axis=1)  # shape: (n_sims_outer, m_steps)

    # Create a single pool and reuse it across all j (avoid respawn each step)
    ctx = get_context("spawn")  # explicit and safe on macOS
    with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as exe:
        for j, t_current in enumerate(time_vector):
            if verbose and (j % 20 == 0 or j == m_steps):
                print(f"  Pricing time step {j}/{m_steps}...")

            m_known_steps = j
            # Build the two 1-D columns we need for this j
            s_col = all_paths[:, j]
            if j > 0:
                sum_known_col = cumsum_S[:, j - 1]
            else:
                sum_known_col = np.zeros(n_sims_outer, dtype=float)

            # Submit by chunks with only 1-D data (much smaller pickles)
            futures = []
            i = 0
            while i < n_sims_outer:
                i_end = min(i + chunk_size, n_sims_outer)
                futures.append(
                    exe.submit(
                        _price_chunk_for_time_light,
                        j, t_current,
                        s_col[i:i_end],
                        sum_known_col[i:i_end],
                        m_known_steps,
                        k, r, sigma, T_maturity, m_steps, omega, n_sims_inner,
                        i, base_seed
                    )
                )
                i = i_end

            # Gather results
            for fut in as_completed(futures):
                i_start, vals = fut.result()
                mtm_matrix[i_start:i_start + len(vals), j] = vals

    return mtm_matrix, time_vector

if __name__ == "__main__":
    s_start = 100
    k = 100
    r = 0.05
    sigma = 0.2
    T_maturity = 1.0
    n_sims_outer = 10_000
    m_steps = 250
    n_sims_inner = 100
    omega = 1

    print("Asian Crude (parallel nested MC)")
    t0 = time.time()

    # Outer paths (vectorized, single-process is fine here)
    all_paths = gbm_path_sim_vectorized(
        s_start, r, sigma, T_maturity, n_sims_outer, m_steps,
        rng=np.random.default_rng(123)
    )

    # Parallel nested pricing
    mtm_matrix_asian, time_vec = calculate_mtm_matrix_asian_parallel(
        all_paths=all_paths, k=k, r=r, sigma=sigma,
        T_maturity=T_maturity, omega=omega, n_sims_inner=n_sims_inner,
        max_workers=None,        # None -> use os.cpu_count()
        chunk_size=1000,         # tune this based on machine
        base_seed=2025,          # reproducible
        verbose=True
    )

    t1 = time.time()
    print(f"Simulation complete. Took {t1 - t0:.2f} seconds.")