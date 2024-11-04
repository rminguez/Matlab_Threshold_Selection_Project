% APPLICATION OF THE HTRESHOLD METHODS TO THE DAILY GLOBAL HISTORICAL CLIMATOLOGY NETWORK (GHCN-DAILY) 
% Version 3.31

% --------------------------------------------------------------------------------
% How to cite:
% 
% Note that the GHCN-Daily dataset itself now has a DOI (Digital Object Identifier)
% so it may be relevant to cite both the methods/overview journal article as well 
% as the specific version of the dataset used.
% 
% The journal article describing GHCN-Daily is:
% Menne, M.J., I. Durre, R.S. Vose, B.E. Gleason, and T.G. Houston, 2012:  An overview 
% of the Global Historical Climatology Network-Daily Database.  Journal of Atmospheric 
% and Oceanic Technology, 29, 897-910, doi:10.1175/JTECH-D-11-00103.1.
% 
% To acknowledge the specific version of the dataset used, please cite:
% Menne, M.J., I. Durre, B. Korzeniewski, S. McNeill, K. Thomas, X. Yin, S. Anthony, R. Ray, 
% R.S. Vose, B.E.Gleason, and T.G. Houston, 2012: Global Historical Climatology Network - 
% Daily (GHCN-Daily), Version 3.31.
% NOAA National Climatic Data Center. http://doi.org/10.7289/V5D21VHZ [2024101719 [yyyymmddhh] UTC].
% --------------------------------------------------------------------------------

% Clear the workspace and close all figures
clear all
close all
clc

% Get the current directory path
actualpath = pwd;
addpath([actualpath '\matlab_functions'])
addpath([actualpath '\data'])

% Read the filelistnames.txt file containing information about files
filelist = readtable('filelistnames.txt', 'Delimiter', ',', 'ReadVariableNames', true);

% Filter files with data availability of 90% or more
filtered_files = filelist(filelist.PercentageAvailable >= 90, :);

% Define significance levels
siglevels = [0.01, 0.05];

% Define the variable names and types for the results table
variableNames = {'Filename', 'Siglevel', 'Langousis_Min1', 'Langousis_Min2', 'Langousis_Min3', 'Langousis_Time', ...
                 'Studentized_Residuals', 'Studentized_Residuals_Time', 'Anderson_Darling', 'Anderson_Darling_Time', ...
                 'CVM', 'CVM_Time'};

% Initialize an empty cell array to store the data
results_cell = {};

% Loop through all filtered files
total_files = height(filtered_files);  % Get the total number of files
for i = 1:total_files
    % Display current progress (file number and total)
    disp(['Processing file ' num2str(i) ' of ' num2str(total_files)]);

    % Get the current filename
    filename = filtered_files.Filename{i};
    disp(['Processing file: ' filename]);
    
    % Check if the data file exists, read it if it does, otherwise throw an error and stop
    if exist(filename, 'file') == 2
        % Read the data file
        data = readtable(filename, 'TreatAsEmpty', {'NA'}, 'Format', '%s%f');
    else
        error('The file %s does not exist. Please check the filename and path.', filename);
    end

    
    % Remove NaN or empty entries
    rowsToDelete = isnan(data{:, 2});
    data(rowsToDelete, :) = [];
    
    % Define precipitation data
    pluviometros.data = data{:,2};

    % Minimum peak distance set to 2 days to consider peaks as independent
    min_peak_distance = 2;
    
    % Extract independent peaks and other variables
    threshold = 0.0;
    n0 = 10;
    [pks_unicos_valid, excedencias_mean_valid, excedencias_weight_valid, pks, ~] = ...
        threshold_peak_extraction(pluviometros.data, threshold, n0, min_peak_distance);
    
    % Apply the optimal threshold detection methods for each significance level
    for siglevel = siglevels
        % 1. Studentized Residuals Method
        tic;  % Start timer for execution time calculation
        [threshold_val_SR, ~, ~, ~] = threshold_studentized_residuals(pks_unicos_valid, ...
            excedencias_mean_valid, excedencias_weight_valid, siglevel);
        time_SR = toc;  % Stop timer for the method
        
        % 2. Langousis Method (local minima)
        tic;  % Start timer for execution time calculation
        threshold_val_MSE = threshold_MSE(pks_unicos_valid, ...
            excedencias_mean_valid, excedencias_weight_valid, n0);
        time_MSE = toc;  % Stop timer for the method
        
        % Store the first three local minima (or NaN if there are less than 3)
        langousis_minima = NaN(1, 3);
        if length(threshold_val_MSE) >= 1
            langousis_minima(1:min(3, length(threshold_val_MSE))) = threshold_val_MSE(1:min(3, length(threshold_val_MSE)));
        end
        
        % 3. Anderson-Darling Method
        tic;  % Start timer for execution time calculation
        threshold_val_AD = threshold_AD(pks, siglevel);
        time_AD = toc;  % Stop timer for the method
        
        % 4. Cramer-Von Mises Method
        tic;  % Start timer for execution time calculation
        threshold_val_CVM = threshold_CVM(pks, siglevel);
        time_CVM = toc;  % Stop timer for the method
         
        % Store the results in the cell array
        results_cell = [results_cell; {filename, siglevel, langousis_minima(1), langousis_minima(2), ...
                                      langousis_minima(3), time_MSE, threshold_val_SR, time_SR, ...
                                      threshold_val_AD, time_AD, threshold_val_CVM, time_CVM}];
        
        % Display the total execution time for the current file and significance level
        disp(['Execution time for ' filename ' with siglevel = ' num2str(siglevel) ': ' num2str(time_MSE+time_SR+time_AD+time_CVM) ' seconds.']);
    end
end

% Convert the cell array to a table
results = cell2table(results_cell, 'VariableNames', variableNames);

% Save the results to a .csv file
writetable(results, 'threshold_results.csv');
disp('Results saved to threshold_results.csv.');
