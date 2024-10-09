# PySurv

PySurv is a Python package for generating and plotting Kaplan-Meier survival curves. It is designed to be simple and intuitive and is inspired from [MatSurv](https://github.com/aebergl/MatSurv).
![KM Curve Example](images/KMCurve.png)
![Summary Example](images/Summary.png)

## Features

- Generate Kaplan-Meier (KM) survival curves
- Calculate hazard ratios (HR) with 95% confidence intervals using Cox PH (default)* or Mantel-Haenszel method
- Option to show or hide confidence intervals (CI) for KM curves
- Flexible plotting with customizable styles, colors, and labels
- Generate synthetic time-to-event data by specifying hazard ratio and censoring rate.

## Installation

You can install PySurv via PyPI:
```
pip install pysurv
```
Make sure you have the following dependencies installed:
```
- numpy
- pandas
- lifelines
- matplotlib
```
These will be installed automatically when you install PySurv using pip.

## Usage

Here’s a detailed example of how to use PySurv for generating data, plotting KM curves, and customizing the output:
```
import pysurv

# Generate synthetic data
data = pysurv.generate_time_to_event_data(n_samples=2000, hazard_ratio=2.0, censoring_rate=0.3)

# Plot Kaplan-Meier curve with custom labels and colors
pysurv.plot_km_curve(
    data, 
    time_col='time', 
    event_col='event', 
    group_col='group', 
    group_labels=('Treatment Group', 'Control Group'),
    title="Survival Analysis of Treatment vs Control",
    y_label="Survival Probability (%)",
    x_label="Time (months)",
    colors=['blue', 'orange'], 
    line_styles=['-', '-'],
    show_ci=True,
    survival_time_point=12,  # Show % survival at 12 months
    return_summary=False
)
```

### Output:

- A Kaplan-Meier curve
- Summary statistics, including median survival time and percentage survival at a given time point
- Hazard ratio and p-value displayed on the plot

## Function Reference

### plot_km_curve

Plots Kaplan-Meier survival curves and calculates hazard ratios, p-values, and confidence intervals.

#### Parameters:
- data: A Pandas DataFrame containing time-to-event data.
- time_col: Column name for the time data.
- event_col: Column name for the event data (1 for event, 0 for censored).
- group_col: Column name for the binary group data.
- group_labels: Labels for the groups (default: ('Group 0', 'Group 1')).
- title: Title for the plot.
- y_label: Label for the y-axis.
- x_label: Label for the x-axis.
- colors: List of colors for the groups (default: ['b', 'r']).
- line_styles: Line styles for the groups (default: ['-', '-']).
- show_ci: Whether to show confidence intervals on the KM curves (default: False).
- method: Method to calculate hazard ratio ('cox'(default)*, 'mantel-haenszel').
- show_inverted_hr: Whether to show inverted hazard ratio (default: False).
- survival_time_point: Time point at which to show percentage survival (default: None).
- return_summary: Whether to return a summary of survival and hazard ratio statistics (default: False).

*Univariable Cox regression model is set to default based on [Frank Harell's advice](https://discourse.datamethods.org/t/when-is-log-rank-preferred-over-univariable-cox-regression/2344). MatSurv uses the Mantel-Haenszel method by default to compute the hazard ratios and log-rank test to compute the p-value.

## ToDo
- Clean up and organize the code
- Support m-ary cases
- Color risk tables below the KM curves

If you’d like to contribute to PySurv, here’s how you can help:

1. Fork the repository.
2. Create a new branch (git checkout -b feature/your-feature).
3. Make your changes and commit them (git commit -am 'Add some feature').
4. Push to the branch (git push origin feature/your-feature).
5. Create a pull request.

Please make sure to update tests as appropriate and follow the style guide used in the project.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Disclaimer
This code is still in the initial stages and may contain bugs. I make no guarantees for the accuracy of this code. By using this code, you agree to take full responsibility for any issues or inaccuracies that may arise.