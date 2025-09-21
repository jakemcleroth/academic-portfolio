import pandas as pd
from sklearn.mixture import GaussianMixture
import numpy as np

# Function to fit GMM and find the best models based on aic, bic and log-likelihood
def fit_gmm(data, n_components=2, n_init=10):
    best_gmm_aic = None
    best_gmm_bic = None
    best_gmm_log_likelihood = None
    best_aic_value = np.inf
    best_bic_value = np.inf
    best_log_likelihood_value = -np.inf

    # Run many times to ensure best parameters are found
    for i in range(n_init):
        # Fit the GMM
        gmm = GaussianMixture(n_components=n_components, covariance_type="diag") # diag for no covariance + can optionally set seed as i here using random_state=i
        gmm.fit(data)
        
        # Don't use if not converged
        if not gmm.converged_:
            print(f"Initialization {i+1}: Did not converge. Skipping.")
            continue
        
        # Metrics
        bic = gmm.bic(data)
        aic = gmm.aic(data)
        log_likelihood = gmm.score(data)
        
        # Update best models
        if aic < best_aic_value:
            best_aic_value = aic
            best_gmm_aic = gmm
        if bic < best_bic_value:
            best_bic_value = bic
            best_gmm_bic = gmm
        if log_likelihood > best_log_likelihood_value:
            best_log_likelihood_value = log_likelihood
            best_gmm_log_likelihood = gmm

    return best_gmm_aic, best_gmm_bic, best_gmm_log_likelihood

# Function to print mixture parameters
def print_model_params(gmm, dataset_name, metric_name):
    if gmm is not None:
        print(f"\nBest Model for {dataset_name} Based on {metric_name}:")
        print(f"Means: {gmm.means_.flatten()}")
        print(f"Covariances: {gmm.covariances_.flatten()}")
        print(f"Weights: {gmm.weights_.flatten()}")
        print(f"BIC: {gmm.bic(ret):.2f}")
        print(f"AIC: {gmm.aic(ret):.2f}")
        print(f"Log-likelihood: {gmm.score(ret):.2f}")
    else:
        print(f"No model converged for {dataset_name} based on {metric_name}.")

if __name__ == "__main__":
    file_path = "/Users/jakemcleroth/Desktop/University/Masters/Modules/Semester 1/VaR/assignment/assignment_1/raw_data/" # can change file path here
    alsi = pd.read_excel(f"{file_path}q5_data_formated.xlsx", sheet_name="ALLSHARE")
    alsi["Return"] = np.log(alsi["JSE ALL SHARE"] / alsi["JSE ALL SHARE"].shift(1))
    alsi_returns = alsi["Return"].dropna().values.reshape(-1, 1)  # Scikit Learn needs a matrix with values and features

    top40 = pd.read_excel(f"{file_path}q5_data_formated.xlsx", sheet_name="TOP40")
    top40["Return"] = np.log(top40["JSE TOP 40"] / top40["JSE TOP 40"].shift(1))
    top40_returns = top40["Return"].dropna().values.reshape(-1, 1)  # Scikit Learn needs a matrix with values and features

    indi25 = pd.read_excel(f"{file_path}q5_data_formated.xlsx", sheet_name="INDUSTRIAL25")
    indi25["Return"] = np.log(indi25["JSE INDUSTRIAL 25"] / indi25["JSE INDUSTRIAL 25"].shift(1))
    indi25_returns = indi25["Return"].dropna().values.reshape(-1, 1)  # Scikit Learn needs a matrix with values and features

    returns = [alsi_returns, top40_returns, indi25_returns]
    names = ["ALSI", "TOP40", "INDI25"]

    # Get best fits for each index
    for ret, name in zip(returns,names):
        print(name)
        best_gmm_aic, best_gmm_bic, best_gmm_log_likelihood = fit_gmm(ret)
        print_model_params(best_gmm_aic, name, "AIC")
        print_model_params(best_gmm_bic, name, "BIC")
        print_model_params(best_gmm_log_likelihood, name, "Log-likelihood")