# xVA: Simulating Exposure for xVA

This project focuses on simulating counterparty credit exposure, with the goal of improving computational efficiency. Different methods are used and compared.

The full write-up is in `xva_project_final.pdf`.

## Notebooks

**Base Simulations**

| File                            | Description                                                                                                                                                                      |
| ------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `simulating_risk_factors.ipynb` | Simulates the underlying risk factors: interest rates via the Hull-White one-factor (HW1F) model calibrated with QuantLib and equity prices via GBM                              |
| `exposure_simulation.ipynb`     | Computes Expected Positive Exposure (EPE) and Potential Future Exposure (PFE) profiles for an interest rate swap, a European option and an Asian option using nested Monte Carlo |

**Improving Performance**

| File                               | Description                                                                                                                         |
| ---------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `antithetics.ipynb`                | Applies antithetic variates to reduce variance in EPE/PFE estimates for European and Asian options, with standard error comparisons |
| `sample_recycling.ipynb`           | Implements the Sample Recycling Method for Asian option exposure                                                                    |
| `vectorization.ipynb`              | Vectorised NumPy implementation of the exposure simulation, replacing Python loops for performance                                  |
| `multiprocessing_asian_options.py` | Parallelised implementation for nested Monte Carlo                                                                                  |
