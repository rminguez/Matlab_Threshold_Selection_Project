function [pks_unicos_valid,excedencias_mean_valid,excedencias_weight_valid, pks, locs] = threshold_peak_extraction(data,threshold,n0,min_peak_distance)
%THRESHOLD_PEAK_EXTRACTION Extracts independent peaks over a threshold from data.
% This function identifies peaks in a dataset that exceed a specified 
% threshold and computes statistics such as mean exceedances, variances, 
% and weights for valid unique peaks. Peaks are considered independent if 
% they are separated by a minimum distance.
%
% INPUTS:
% data               - Input time series or data array
% threshold          - Threshold above which peaks are extracted
% n0                 - Minimum number of exceedances required for valid computation
% min_peak_distance  - Minimum distance between two peaks (in data points)
%
% OUTPUTS:
% pks_unicos_valid           - Valid unique peaks after removing NaN values
% excedencias_mean_valid     - Mean exceedances for valid peaks
% excedencias_weight_valid   - Weights based on exceedance variance for valid peaks
% pks                       - All detected peaks
% locs                      - Indices of the detected peaks in the data
%
% This function implements the Peaks Over Threshold (POT) method with the 
% requirement of independent peaks to avoid biases caused by dependent extremes.

    % Find peaks that exceed the threshold with a specified minimum distance
    [pks, locs] = findpeaks(max(data-threshold,0), ...
                        'MinPeakDistance', min_peak_distance);
    
    %%%%%%%%%%%%%%%%%%

    % Eliminate duplicate peaks and sort them in ascending order
    pks_unicos = unique(pks);

    % Initialize arrays to store mean exceedances, variances, and weights
    excedencias_mean = zeros(length(pks_unicos), 1);
    excedencias_var = zeros(length(pks_unicos), 1);
    excedencias_weight = zeros(length(pks_unicos), 1);

    % Loop through each unique peak and calculate mean exceedances, variances, and weights
    for i = 1:length(pks_unicos)
        % Define the current unique peak
        pico_actual = pks_unicos(i);

        % Calculate the exceedances for peaks greater than the current unique peak
        excedencias = pks(pks > pico_actual);

        % If there are enough exceedances (greater than or equal to n0)
        if length(excedencias) >= n0
            % Compute the mean exceedance (adjusted by the current peak)
            excedencias_mean(i) = mean(excedencias) - pico_actual;
            % Compute the variance of the exceedances
            excedencias_var(i) = var(excedencias);
            % Compute the weight as the number of exceedances divided by the variance
            excedencias_weight(i) = length(excedencias) / excedencias_var(i);
        else
            % If fewer than n0 exceedances, truncate the arrays and stop the loop
            excedencias_mean = excedencias_mean(1:(i-1));
            excedencias_var = excedencias_var(1:(i-1));
            excedencias_weight = excedencias_weight(1:(i-1));
            break
        end
    end

    % Trim the list of unique peaks to match the number of valid exceedances
    pks_unicos = pks_unicos(1:length(excedencias_weight));

    % Remove any NaN values from the peak and exceedance data to avoid issues in regression
    valid_indices = ~isnan(pks_unicos) & ~isnan(excedencias_mean) & ~isnan(excedencias_weight);
    pks_unicos_valid = pks_unicos(valid_indices);
    excedencias_mean_valid = excedencias_mean(valid_indices);
    excedencias_weight_valid = excedencias_weight(valid_indices);

end

