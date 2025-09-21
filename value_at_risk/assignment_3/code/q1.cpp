#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

using namespace std;

// TO RUN CODE: (It works nicely in VSCODE if you have it)
// C++ compiler is already installed on Mac (not Windows)
// Step 1. Type: "g++ -std=c++11 q1.cpp -o q1" into command line/terminal to compile code
// Step 2. Type: "./q1" to run code
// VaR Values will output to terminal

// Function that prints the value at risk to command line. Takes in n for number of sims, h for days, a vector levels for alph levels and then an array norm_sims for the standard normal simulations
void value_at_risk_cv(double vol, int n, int h, const vector<vector<double>> &norm_sims, const vector<double> &levels) {

    double daily_vol = vol / sqrt(250.0);

    // Calculate daily returns from normal sims
    vector<vector<double>> returns_cv(n, vector<double>(h));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < h; ++j)
            returns_cv[i][j] = norm_sims[i][j] * daily_vol; // Just multiplying by the vol

    // Calculate h-day returns i.e. the sum of each row
    vector<double> returns_cv_hday(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < h; ++j)
            sum += returns_cv[i][j];
        returns_cv_hday[i] = sum;
    }

    // Sort the returns (so that we can compute the quantile easily)
    sort(returns_cv_hday.begin(), returns_cv_hday.end());

    // Loop to calc VaR for each alpha i.e each level in levels
    cout << "Value at Risk (Constant Vol):" << endl;
    for (double level : levels) {
        size_t idx = static_cast<size_t>(round(level * n)); // Position of quantile in simulated data
        if (idx >= n) idx = n - 1; // Rounding might cause idx = n
        double var = -returns_cv_hday[idx]; // Negative alpha quantile
        cout << "VaR at " << level * 100 << "%: " << var << endl;
    }
}

// Function that print value at risk for agarch to command line. Takes in agarch parameters, first return, n for number of sims, h for days, a vector levels for alph levels and then an array norm_sims for the standard normal simulations
void value_at_risk_agarch(double omega, double lambda, double alpha, double beta, double current_ret, int n, int h, const vector<vector<double>> &norm_sims, const vector<double> &levels) {
    // Parameters
    double vol_unc_daily = sqrt((omega + (pow(lambda,2))*alpha)/(1-(alpha+beta)));

    // Initialize vol and return matrices
    vector<vector<double>> vol_agarch(n, vector<double>(h));
    vector<vector<double>> returns_agarch(n, vector<double>(h));

    // Set initial values (t = 0)
    for (int i = 0; i < n; ++i) {
        vol_agarch[i][0] = sqrt(omega + alpha * pow(current_ret - lambda, 2) + beta * pow(vol_unc_daily, 2));
        returns_agarch[i][0] = vol_agarch[i][0] * norm_sims[i][0];
    }

    // A-GARCH Simulation for returns and vol
    for (int t = 1; t < h; ++t) {
        for (int i = 0; i < n; ++i) {
            vol_agarch[i][t] = sqrt(omega +
                                    alpha * pow(returns_agarch[i][t - 1] - lambda, 2) +
                                    beta * pow(vol_agarch[i][t - 1], 2));
            returns_agarch[i][t] = vol_agarch[i][t] * norm_sims[i][t];
        }
    }

    // h-day A-GARCH returns
    vector<double> returns_agarch_hday(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int t = 0; t < h; ++t)
            sum += returns_agarch[i][t];
        returns_agarch_hday[i] = sum;
    }

    // Sort the returns
    sort(returns_agarch_hday.begin(), returns_agarch_hday.end());

    // Loop to calc VaR for each alpha i.e each level in levels
    cout << "Value at Risk (A-GARCH):" << endl;
    for (double level : levels) {
        size_t idx = static_cast<size_t>(round(level * n)); // Position of quantile in simulated data. size_t is specifically used for sizes/indices
        if (idx >= n) idx = n - 1; // Rounding might cause idx = n
        double var = -returns_agarch_hday[idx]; // Negative alpha quantile
        cout << "VaR at " << level * 100 << "%: " << var << endl;
    }
}

int main() {
    // GLOBAL VARIABLES
    int n = 1000000; //no. of sims
    int h = 5; //no. of days
    vector<double> levels = {0.001, 0.01, 0.05, 0.1}; // alpha levels

    // Set up for simulation
    random_device rd; // Random device object - non-deterministic seed
    mt19937 gen(rd()); // Mersenne Twister
    normal_distribution<> d(0.0, 1.0); // Normal distribution object that uses gen to simulate normal observations

    // Simulate standard normal values in an n x h structure
    vector<vector<double>> norm_sims(n, vector<double>(h));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < h; ++j)
            norm_sims[i][j] = d(gen);


    // Constant Volatility VaR
    double vol = 0.25; // vol for constant volatility VaR
    value_at_risk_cv(vol, n, h, norm_sims, levels); // Call function to output VaR


    // A-GARCH VaR
    double omega = 4*pow(10,-6);
    double alpha = 0.06;
    double lambda = 0.01;
    double beta = 0.9;
    double current_ret_pos = 0.1;
    double current_ret_neg = -0.1;
    value_at_risk_agarch(omega, lambda, alpha, beta, current_ret_pos,n, h, norm_sims, levels);
    value_at_risk_agarch(omega, lambda, alpha, beta, current_ret_neg,n, h, norm_sims, levels);
    return 0;
}