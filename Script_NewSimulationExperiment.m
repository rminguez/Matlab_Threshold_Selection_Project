% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT to evaluate threshold selection methods for different parameter values
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

% Define paths
actualpath = pwd;
addpath([actualpath '\matlab_functions']);

warning('off', 'all'); % Suppress warnings
rng(42); % Set random seed for reproducibility

% Experiment parameters
cdfumbralY0_vals = [0.99, 0.995, 0.999];
k_vals = [-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15];
sigma_vals = [0.5, 1, 2];
siglevel_vals = [0.01, 0.05];
rond_vals = [1, 2, 3]; % Rounding precision options
n_years = 40; % Simulation length in years
n = 365.25 * n_years; % Total number of days
num_runs = 10; % Number of runs per parameter combination

% Calculate the total number of iterations
total_iterations = length(cdfumbralY0_vals) * length(k_vals) * length(sigma_vals) * ...
                   length(siglevel_vals) * length(rond_vals) * num_runs;
iteration_counter = 0; % Initialize the iteration counter

% Initialize results storage
results = [];
varnames = {'Iteration', 'Run', 'CDF_Threshold', 'k', 'sigma', 'Rounding', 'Real_Threshold', 'SigLevel', ...
            'Studentized_Residuals', 'SR_Time', ...
            'Langousis_Min1', 'Langousis_Min2', 'Langousis_Min3', 'Langousis_Time', ...
            'Anderson_Darling', 'AD_Time', ...
            'Cramer_Von_Mises', 'CVM_Time'};

% Loop over parameter combinations
for cdfumbralY0 = cdfumbralY0_vals
    for k = k_vals
        for sigma = sigma_vals
            for run = 1:num_runs
                % Simulate daily data
                ySim = max(normrnd(0, 1, n, 1), 0); % Gaussian with zero truncation
                threshold0 = norminv(cdfumbralY0); % Real threshold for Gaussian data
                
                % Adjust tail with Generalized Pareto Distribution (GPD)
                cola = (ySim >= threshold0);
                if sum(cola) > 0
                    zcola = ySim(cola);
                    probscaled = (normcdf(zcola) - cdfumbralY0) / (1 - cdfumbralY0);
                    ySim(cola) = gpinv(probscaled, k, sigma, threshold0); % Apply GPD
                end
                
                % Loop over rounding options
                for rond = rond_vals
                    pluviometros.data = round(ySim, rond); % Round to match daily precision
                    
                    % Extract peaks
                    threshold = 0.0;
                    n0 = 10;
                    min_peak_distance = 1;
                    [pks_unicos_valid, excedencias_mean_valid, excedencias_weight_valid, pks] = ...
                        threshold_peak_extraction(pluviometros.data, threshold, n0, min_peak_distance);
                    
                    % Loop over significance levels
                    for siglevel = siglevel_vals
                        % Increment the iteration counter
                        iteration_counter = iteration_counter + 1;
                        
                        % Display progress
                        fprintf('Iteration %d/%d - Run: %d, CDF: %.3f, k: %.2f, sigma: %.1f, Rounding: %d, SigLevel: %.2f\n', ...
                                iteration_counter, total_iterations, run, cdfumbralY0, k, sigma, rond, siglevel);
                        
                        % 1. Studentized Residuals Method
                        tic;
                        [threshold_val_SR, ~, ~, ~] = threshold_studentized_residuals(pks_unicos_valid, ...
                            excedencias_mean_valid, excedencias_weight_valid, siglevel);
                        time_SR = toc;
                        
                        % 2. Langousis Method
                        tic;
                        threshold_val_MSE = threshold_MSE(pks_unicos_valid, ...
                            excedencias_mean_valid, excedencias_weight_valid, n0,0.999);
                        time_MSE = toc;
                        
                        % Store the first three local minima (or NaN if fewer)
                        langousis_minima = NaN(1, 3);
                        if length(threshold_val_MSE) >= 1
                            langousis_minima(1:min(3, length(threshold_val_MSE))) = threshold_val_MSE(1:min(3, length(threshold_val_MSE)));
                        end
                        
                        % 3. Anderson-Darling Method
                        tic;
                        threshold_val_AD = threshold_AD(pluviometros.data, siglevel);
                        time_AD = toc;
                        
                        % 4. Cramer-Von Mises Method
                        tic;
                        threshold_val_CVM = threshold_CVM(pluviometros.data, siglevel);
                        time_CVM = toc;
                        
                        % Append results
                        results = [results; {iteration_counter, run, cdfumbralY0, k, sigma, rond, threshold0, siglevel, ...
                                             threshold_val_SR, time_SR, ...
                                             langousis_minima(1), langousis_minima(2), langousis_minima(3), time_MSE, ...
                                             threshold_val_AD, time_AD, ...
                                             threshold_val_CVM, time_CVM}];
                        
                    end
                end
            end
        end
    end
end

% Save results to CSV
results_table = cell2table(results, 'VariableNames', varnames);
writetable(results_table, 'simulation_results_with_rounding_and_runs.csv');
disp('All results saved to simulation_results_with_rounding_and_runs.csv');

warning('on', 'all'); % Re-enable warnings
