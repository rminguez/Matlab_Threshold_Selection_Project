# Matlab Threshold Selection Project

This repository contains Matlab scripts and functions for threshold selection in precipitation data, specifically for automating threshold selection using various statistical tests. These methods are applied to precipitation data extracted from the **Global Historical Climatology Network - Daily (GHCN-Daily)** dataset, Version 3.31.

## Dataset Citation
To acknowledge the specific dataset version used, please cite as follows:
- Menne, M.J., I. Durre, R.S. Vose, B.E. Gleason, and T.G. Houston, 2012: An overview of the Global Historical Climatology Network-Daily Database. *Journal of Atmospheric and Oceanic Technology*, 29, 897-910, doi:10.1175/JTECH-D-11-00103.1.
- Menne, M.J., et al., 2012: Global Historical Climatology Network - Daily (GHCN-Daily), Version 3.31. *NOAA National Climatic Data Center*. [https://doi.org/10.7289/V5D21VHZ](https://doi.org/10.7289/V5D21VHZ).

## Project Structure
- **data**: Contains CSV files and any data extracted from the GHCN-Daily dataset.
- **figures**: All generated plots, including histograms, boxplots, and threshold comparison visualizations.
- **matlab_functions**: Core Matlab functions used for statistical analysis and threshold selection.

## Scripts

### `Script_Download_Rainfall_Database.m`
This script downloads the GHCN-Daily dataset's metadata and filters stations with over 110 years of data. It processes `.dly` files and extracts daily precipitation records for each station, storing them in `data`.

### `Script_GHCN_DAILY_Database.m`
This script automates threshold selection on the GHCN-Daily dataset, applying several statistical methods:
- **Studentized Residuals**
- **Langousis Method**
- **Anderson-Darling Test**
- **Cramer-Von Mises Test**

For each file with over 90% data availability, the script iterates through significance levels of 0.01 and 0.05 to determine optimal thresholds based on each method.

### `Script_GHCN_DAILY_AnalysesResults.m`
Analyzes and visualizes the results from the threshold selection process. Generates summary statistics, boxplots, and performs additional statistical tests to compare the methods. Key plots include:
- **Threshold Comparison Across Methods**: A boxplot of selected thresholds.
- **CPU Time Comparison**: A boxplot comparing computation times for each method.

### `Script_New_PRCP_ASN00021043.m`
Illustrates the application of these methods on a single station (`PRCP_ASN00021043`). It generates individual threshold selections and visualizes the precipitation time series with independent peaks marked.

## Core Functions

### `threshold_peak_extraction.m`
Identifies independent peaks over a given threshold, computes mean exceedances, variances, and weights for valid unique peaks. This function is the basis for all threshold selection methods by ensuring independent data points.

### `threshold_studentized_residuals.m`
Selects an optimal threshold using the Chi-squared and studentized residuals tests, aiming for independence in the extreme values used in analysis.

### `threshold_MSE.m`
Implements the **Langousis Method**, calculating candidate thresholds based on local minima of a spline fit to the data.

### `threshold_AD.m`
Applies the **Anderson-Darling Test** to select the optimal threshold by fitting a Generalized Pareto Distribution to the data's excesses.

### `threshold_CVM.m`
Uses the **Cramer-Von Mises Test** to determine the best threshold for data exceeding unique values, fitting a Generalized Pareto Distribution and selecting based on CVM statistics.

### `RWLSfit.m`
Performs robust weighted least squares regression, used for calculating residuals in the **Studentized Residuals Method**.

## Compatibility

This project has been tested using MATLAB 2016. Compatibility with other versions of MATLAB has not been fully verified, so we recommend using MATLAB 2016 for best results.

## Usage Instructions

1. **Download and Install Matlab**: Ensure Matlab is installed and configured on your system.
2. **Clone this Repository**:
   ```bash
   git clone https://github.com/rminguez/Matlab_Threshold_Selection_Project.git
