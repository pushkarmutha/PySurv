import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
import matplotlib.pyplot as plt
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

# Function to plot KM curves
def plot_km_curve(data, time_col='time', event_col='event', group_col='group', 
                  group_labels=None, title=None, 
                  y_label='Survival Probability', x_label='Time (months)', colors=['r', 'b'], line_styles=None, fontsize=18, linewidth=2.5,
                  show_ci=False, show_inverted_hr=False, survival_time_points=None, return_summary=False, savepath=None, **kwargs):
    """
    Plots Kaplan-Meier survival curves and displays hazard ratio, p-value, and confidence intervals.
    
    Parameters:
    data (pd.DataFrame): The time-to-event dataset.
    time_col (str): Column name for time data.
    event_col (str): Column name for event data (1 for event occurred, 0 for censored).
    group_col (str): Column name for group data.
    group_labels (tuple): Labels for the two groups (default: None).
    title (str): Title for the plot (default: None).
    y_label (str): Label for the y-axis (default: 'Survival Probability').
    x_label (str): Label for the x-axis (default: 'Time (months)').
    colors (list): List of colors to use for the groups (default: ['r', 'b']). If more than two groups, please manually provide a list of colors.
    line_styles (list): List of line styles to use for the groups (default: '-').
    fontsize (int): Font size for text on KM curve including title, axis labels, and risk tables (default: 18).
    linewidth (float): Line width of KM curves (default: 2.5).
    show_ci (bool): Whether to show confidence intervals on KM curves (default: False).
    show_inverted_hr (bool): Whether to show inverted hazard ratio (default: False).
    survival_time_points (list): One or more time point(s) at which to estimate percentage survival (default: None).
    return_summary (bool): Whether to return a summary of survival and hazard ratio statistics (default: False).
    savepath (str): Complete path (including filename and extension) to save the KM curve plot (default: None). 
    **kwargs: Additional matplotlib arguments to pass for plotting KM curves.
    Returns:
    survival_summary: If return_summary=True, Pandas dataframe with median survival and % patients alive at specified timepoint
    hr_summary: If return_summary=True, Pandas dataframe with hazard ratio, confidence interval, p-value and test statistic
    """
    if  data[group_col].nunique()!=2:
        print('Please explicitly provide a list of "colors"')
        
    groups = sorted(data[group_col].unique())
    
    if line_styles is None:
        line_styles = ['-'] * len(groups)

    if group_labels is None:
        group_labels = groups

    plt.figure(figsize=(12, 8), facecolor='white')
    
    ax = plt.subplot(111)
    survival_percentages = []
    kmfs = []

    for i, group in enumerate(groups):
        kmf = KaplanMeierFitter()
        group_data = data[data[group_col] == group]
        kmf.fit(group_data[time_col], event_observed=group_data[event_col], label=group_labels[i])
        kmf.plot_survival_function(show_censors=True, censor_styles={"marker": "|", "ms":10, "mew": 1.25}, ci_show=show_ci, ci_alpha=0.15, color=colors[i], linestyle=line_styles[i], ax=ax, fontsize=fontsize, linewidth=linewidth, **kwargs)
        kmfs.append(kmf)

        # Record median survival and percentage survival at a specific time point if provided
        median_survival = kmf.median_survival_time_
        if survival_time_points is None:
            survival_percentages.append((group_labels[i], median_survival))
        else:
            survival_at_tps = np.array(kmf.predict(survival_time_points))*100
            if len(survival_time_points)==1:
                survival_percentages.append((group_labels[i], median_survival, survival_at_tps))
            else:
                survival_percentages.append((group_labels[i], median_survival, *survival_at_tps))

    hr = None
    ci_lower = None
    ci_upper = None
    p_value = None
    
    # Fit Cox Proportional Hazards model to calculate hazard ratio and p-value
    cph = CoxPHFitter()
    cph.fit(data[[group_col, time_col, event_col]], duration_col=time_col, event_col=event_col)
    hr = cph.hazard_ratios_[group_col]
    ci_lower, ci_upper = np.exp(cph.confidence_intervals_.loc[group_col])
    p_value = cph.summary.loc[group_col, 'p']
    test_statistic = cph.summary.loc[group_col, 'z']

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
        plt.text(0.0125, 0.025, f"HR: {hr:.2f} ({ci_lower:.2f}-{ci_upper:.2f})\np-value: {p_value}",
                 horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes, 
                 bbox=dict(facecolor='white', alpha=0.3, edgecolor='none'), fontsize=fontsize-2)
    remove_spines(ax,['top', 'right'])
    plt.xlabel(x_label, fontsize=fontsize)
    plt.ylabel(y_label, fontsize=fontsize)
    if title:
        plt.title(r'$\bf{'+str(title).replace(' ', r'\ ')+'}$', fontsize=fontsize)
    plt.grid(False)

    plt.legend(fontsize=fontsize)
    add_at_risk_counts(*kmfs, ax=ax, fontsize=fontsize)
    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout()
    if savepath is not None:
        plt.savefig(savepath, facecolor='white')
    plt.show()
    
    # Print summary table
    print("\nSummary Table:")
    if survival_time_points is None:
        survival_summary = pd.DataFrame(survival_percentages, columns=['Group', 'Median Survival Time'])
    else:
        surv_tps_cols = [f'% survival at {i}' for i in survival_time_points]
        column_names = ['Group', 'Median Survival Time', *surv_tps_cols]
        survival_summary = pd.DataFrame(survival_percentages, columns=column_names)
    print(survival_summary.round(2).to_string())

    # Print Hazard Ratio Summary
    if hr is not None:
        print("\nHazard Ratio Summary:")
        print(f"Hazard ratio computed using Cox univariable regression on {group_col} variable")
        hr_summary = pd.DataFrame({
            'Hazard Ratio': [hr],
            '95% CI Lower': [ci_lower],
            '95% CI Upper': [ci_upper],
            'P-value': [p_value_exact],
            'Test Statistic': [test_statistic]
        })
        print(hr_summary.round({"Hazard Ratio":2, "95% CI Lower":2, "95% CI Upper":2, "Test Statistic":2}).to_string())
        print("\n\nUnder the proportional hazards assumption, Cox regression is equivalent to log-rank test. See https://www.fharrell.com/post/logrank/ for more info.")

    if return_summary:
        return survival_summary, hr_summary