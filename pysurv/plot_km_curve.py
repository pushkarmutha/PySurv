import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
import matplotlib.pyplot as plt
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts, remove_spines

# Function to generate random time-to-event data with censoring
def generate_time_to_event_data(n_samples=200, hazard_ratio=2.0, censoring_rate=0.3, seed=42):
    np.random.seed(seed)
    
    # Group assignment
    group = np.random.choice([0, 1], size=n_samples)
    
    # Generate survival times using an exponential distribution
    # The hazard ratio is reflected by changing the scale parameter
    baseline_hazard = 0.05
    survival_times = np.where(
        group == 0,
        np.random.exponential(1 / baseline_hazard, n_samples),
        np.random.exponential(1 / (baseline_hazard * hazard_ratio), n_samples)
    )
    
    # Generate random censoring times
    censoring_times = np.random.exponential(1 / (baseline_hazard * censoring_rate), n_samples)
    
    # Determine observed times and censoring status
    observed_times = np.minimum(survival_times, censoring_times)
    event_observed = survival_times <= censoring_times
    
    # Create a DataFrame with the generated data
    data = pd.DataFrame({
        'time': observed_times,
        'event': event_observed,
        'group': group
    })
    
    return data

# Function to perform log-rank and Mantel-Haenszel test
def mantel_haenszel_test(data, time_col='time', event_col='event', group_col='group'):
    groups = data[group_col].unique()
    assert len(groups) == 2, "This implementation only supports two groups"
    
    # Extract data for each group
    group_0_data = data[data[group_col] == groups[0]].sort_values(by=time_col)
    group_1_data = data[data[group_col] == groups[1]].sort_values(by=time_col)
    
    # Get unique event times from both groups
    unique_times = np.unique(data[time_col])
    n_times = len(unique_times)
    
    # Allocate matrices for calculations
    mf = np.zeros((n_times, 2))  # Observed failures
    nf = np.zeros((n_times, 2))  # Number at risk
    ef = np.zeros((n_times, 2))  # Expected number of failures
    
    # Calculate number at risk and observed failures for each group
    for i, t in enumerate(unique_times):
        nf[i, 0] = np.sum(group_0_data[time_col] >= t)
        nf[i, 1] = np.sum(group_1_data[time_col] >= t)
        mf[i, 0] = np.sum((group_0_data[time_col] == t) & (group_0_data[event_col]))
        mf[i, 1] = np.sum((group_1_data[time_col] == t) & (group_1_data[event_col]))

    # Calculate expected number of failures
    nf_sum = np.sum(nf, axis=1)
    mf_sum = np.sum(mf, axis=1)
    for i in range(2):
        ef[:, i] = (nf[:, i] / nf_sum) * mf_sum
    
    results = logrank_test(group_0_data['time'], group_1_data['time'], group_0_data['event'], group_1_data['event'])
    chi2_stat, p_value = results.test_statistic, results.p_value
    # Calculate hazard ratio and 95% confidence interval for Mantel-Haenszel
    hr_mh = (np.sum(mf[:, 1]) / np.sum(ef[:, 1])) / (np.sum(mf[:, 0]) / np.sum(ef[:, 0]))
    log_hr_se = np.sqrt(1 / np.sum(ef[:, 0]) + 1 / np.sum(ef[:, 1]))
    ci_lower = np.exp(np.log(hr_mh) - 1.96 * log_hr_se)
    ci_upper = np.exp(np.log(hr_mh) + 1.96 * log_hr_se)
    
    return chi2_stat, p_value, hr_mh, ci_lower, ci_upper

# Function to plot KM curves
def plot_km_curve(data, time_col='time', event_col='event', group_col='group', 
                  group_labels=('Group 0', 'Group 1'), title=None, 
                  y_label='Survival Probability', x_label='Time (months)', colors=None, line_styles=None, 
                  show_ci=False, method='cox', show_inverted_hr=False, survival_time_point=None, return_summary=False):
    """
    Plots Kaplan-Meier survival curves and displays hazard ratio, p-value, and confidence intervals.
    
    Parameters:
    data (pd.DataFrame): The time-to-event dataset.
    time_col (str): Column name for time data.
    event_col (str): Column name for event data (1 for event occurred, 0 for censored).
    group_col (str): Column name for group data.
    group_labels (tuple): Labels for the two groups.
    title (str): Title for the plot.
    y_label (str): Label for the y-axis.
    x_label (str): Label for the x-axis.
    colors (list): List of colors to use for the groups.
    line_styles (list): List of line styles to use for the groups.
    show_ci (bool): Whether to show confidence intervals on KM curves.
    method (str): Method for calculating hazard ratio ('cox'(default), 'mantel-haenszel').
    show_inverted_hr (bool): Whether to show inverted hazard ratio.
    survival_time_point (float): Time point at which to show percentage survival.
    
    Returns:
    summary_table: Pandas dataframe with median survival and % patients alive at specified timepoint
    hr_summary: Pandas dataframe with hazard ratio, confidence interval, p-value and test statistic
    """
    
    if colors is None:
        colors = ['b', 'r']
    if line_styles is None:
        line_styles = ['-', '-']
        
    groups = data[group_col].unique()
    plt.figure(figsize=(12, 8))
    
    ax = plt.subplot(111)
    survival_percentages = []
    kmfs = []

    for i, group in enumerate(groups):
        kmf = KaplanMeierFitter()
        group_data = data[data[group_col] == group]
        kmf.fit(group_data[time_col], event_observed=group_data[event_col], label=group_labels[i])
        kmf.plot_survival_function(show_censors=True, censor_styles={"marker": "|", "ms":6}, ci_show=show_ci, ci_alpha=0.15, color=colors[i], linestyle=line_styles[i], ax=ax)
        kmfs.append(kmf)
        # Add censoring marks to the curve as dashes
        #censor_times = group_data[time_col][(~group_data[event_col])]
        #plt.vlines(censor_times, ymin=[kmf.predict(t) - 0.01 for t in censor_times], ymax=[kmf.predict(t) + 0.01 for t in censor_times], color=colors[i], linestyles='solid', alpha=0.6)

        # Record median survival and percentage survival at a specific time point if provided
        median_survival = kmf.median_survival_time_
        survival_percentages.append((group_labels[i], median_survival, kmf.predict(survival_time_point)*100 if survival_time_point else None))

    hr = None
    ci_lower = None
    ci_upper = None
    p_value = None
    
    if method == 'cox':
        # Fit Cox Proportional Hazards model to calculate hazard ratio and p-value
        cph = CoxPHFitter()
        cph.fit(data[[group_col, time_col, event_col]], duration_col=time_col, event_col=event_col)
        hr = cph.hazard_ratios_[group_col]
        ci_lower, ci_upper = np.exp(cph.confidence_intervals_.loc[group_col])
        p_value = cph.summary.loc[group_col, 'p']
        test_statistic = cph.summary.loc[group_col, 'z']
    elif method == 'mantel-haenszel':
        # Perform log-rank or Mantel-Haenszel test
        test_statistic, p_value, hr, ci_lower, ci_upper = mantel_haenszel_test(data, time_col, event_col, group_col)

    if show_inverted_hr and hr is not None:
        hr = 1 / hr
        ci_lower_inv = 1 / ci_upper
        ci_upper_inv = 1 / ci_lower
        ci_lower = ci_lower_inv
        ci_upper = ci_upper_inv
    
    if p_value < 0.0001:
        p_value_exact = p_value.copy()
        p_value = "p < 0.0001"
    else:
        p_value_exact = p_value.copy()
        p_value = str(round(p_value, 4))

    # Display hazard ratio, confidence interval, and p-value inside the plot near bottom left
    if hr is not None and ci_lower is not None and ci_upper is not None and p_value is not None:
        plt.text(0.05, 0.05, f"HR: {hr:.2f} ({ci_lower:.2f}-{ci_upper:.2f})\np-value: {p_value}",
                 horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes, fontsize=12, 
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
    remove_spines(ax,['top', 'right'])
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if title:
        plt.title(title)
    plt.grid(False)

    plt.legend()
    add_at_risk_counts(*kmfs, ax=ax)
    plt.subplots_adjust(bottom=0.3)
    plt.show()
    
    # Print summary table
    print("\nSummary Table:")
    summary_table = pd.DataFrame(survival_percentages, columns=['Group', 'Median Survival Time', f'% Alive at {survival_time_point}'])
    print(summary_table)

    # Print Hazard Ratio Summary
    if hr is not None:
        print("\nHazard Ratio Summary:")
        print(f"Hazard ratio computed using {method} method")
        hr_summary = pd.DataFrame({
            'Hazard Ratio': [hr],
            '95% CI Lower': [ci_lower],
            '95% CI Upper': [ci_upper],
            'P-value': [p_value_exact],
            'Test Statistic': [test_statistic]
        })
        print(hr_summary)
    if return_summary:
        return summary_table, hr_summary