"""
The following python script calculates the returns after fees as well as fees earned for a given hedge fund. 
The output is saved in excel spreadsheets.
Assumption: For this task, we assume that the returns given are always before fees.
"""

import pandas as pd


# Annual Fee = Management Fee + {Max[0, Incentive Fee × (Gross Return above HWM − Management Fee − Hurdle Rate)]}
def fund_fees_calculator(capital, returns, management, perf_fee=0.0, hurdle=0.0, has_hwm=False):

    """
    Calculate the returns and rand amount after fees for a given hedge fund

    Args:
        capital: amount invested (set to R1 000 000 but can be changed)
        returns: csv of monthly returns
        perf_fee: performance fee of hedge fund
        hurdle: hurdle rate of hedge fund
        has_hwm: true if hedge fund has a high water mark, false otherwise

    Returns:
        pd.DataFrame: columns for year, starting capital, management fees (R), ending capital (before and after management fees), annual return (before and after management fees), performance fee (R), 
        ending capital after all fees, net return after all fees and total fees earned
    """

    management_fee_rate = management / 12  # Monthly
    capital_tracker = capital
    hwm = capital  # Initial High Water Mark
    results = []

    for year in returns.index:
        starting_cap = capital_tracker # starting capital for the year
        gross_return = 1.0

        # Apply monthly returns and management fee
        management_fee = 0.0
        for month in returns.columns:
            monthly_ret = returns.loc[year, month]
            if pd.isna(monthly_ret):
                continue

            capital_tracker *= (1 + monthly_ret) # increase capital by return amount each month
            gross_return *= (1 + monthly_ret) # keep track of returns without fees
            management_fee = management_fee + capital_tracker*management_fee_rate # management fees earned in rands
            capital_tracker *= (1 - management_fee_rate) # decrease capital by management fee after return for the month has been added
            
        # Annual gross return
        annual_gross_return = gross_return - 1

        # Save capital value amount without performance fees for later calculationg of returns without performance fee but with management fees
        capital_tracker_without_performance_fee = capital_tracker
        
        # Compute excess return for performance fee
        performance_fee = 0.0
        if perf_fee > 0:
            # Calculate return over HWM
            excess_return_capital_amount = capital_tracker - hwm - hurdle*starting_cap # capital_tracker has already had management fees subtracted, hence its not included here
            # Apply the performance fees
            if excess_return_capital_amount > 0:
                performance_fee = perf_fee * excess_return_capital_amount
                capital_tracker = capital_tracker - performance_fee

        # Update HWM
        if capital_tracker > hwm:
            hwm = capital_tracker

        results.append({
            "Year": year,
            "Starting Capital": starting_cap,
            "Management Fees (R)": management_fee,
            "Ending Capital (pre-management fee)": starting_cap*(1+annual_gross_return),
            "Ending Capital (post-management fee)": capital_tracker_without_performance_fee,
            "Annual Return (pre-management fee)": annual_gross_return,
            "Annual Return (post-management fee)": capital_tracker_without_performance_fee/starting_cap - 1,
            "Performance Fees (R)": performance_fee,
            "Ending Capital (after all fees)": capital_tracker,
            "Net Return (after all fees)": (capital_tracker - starting_cap) / starting_cap,
            "Total Fees Earned (R)": performance_fee + management_fee
        })

    return pd.DataFrame(results).set_index("Year")

def main():
    # Get data
    matrix_returns = pd.read_csv("data/matrix.csv", index_col='Year')
    matrix_returns = matrix_returns.sort_index(ascending=True)/100 # matrix values weren't in decimals
    keystone_returns = pd.read_csv("data/keystone.csv", index_col='Year')
    polar_returns = pd.read_csv("data/polar.csv", index_col='Year')
    rocksolid_returns = pd.read_csv("data/rocksolid.csv", index_col='Year')
   
    print("Matrix: ")
    matrix = fund_fees_calculator(
        capital=1000000,
        returns=matrix_returns,
        management=0.0125,
        perf_fee=0.2,
        hurdle=0.0,
        has_hwm=True
    )
    matrix.to_excel("results/matrix_results.xlsx")
    print(matrix)

    print("Polar: ")
    polar = fund_fees_calculator(
        capital=1000000,
        returns=polar_returns,
        management=0.026,
        perf_fee=0.2,
        hurdle=0.07,
        has_hwm=True
    )
    polar.to_excel("results/polar_results.xlsx")
    print(polar)

    print("Keystone: ")
    keystone = fund_fees_calculator(
        capital=1000000,
        returns=keystone_returns,
        management=0.0105,
        perf_fee=0.15,
        hurdle=0.07,
        has_hwm=True
    )
    keystone.to_excel("results/keystone_results.xlsx")
    print(keystone)

    print("Rock Solid: ")
    rocksolid = fund_fees_calculator(
        capital=1000000,
        returns=rocksolid_returns,
        management=0.0125,
        perf_fee=0.10,
        hurdle=0.07,
        has_hwm=True
    )
    rocksolid.to_excel("results/rocksolid_results.xlsx")
    print(rocksolid)



if __name__ == "__main__":
    main()
    
