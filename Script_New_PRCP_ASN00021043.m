%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %   SCRIPT to optimize the selection of the threshold for Record PRCP_ASN00021043
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
actualpath=pwd;

% Add folders to the paths
addpath([actualpath '\matlab_functions'])

%% Plotting the data
% 
graficos = 1; % Enable/disable graphics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXTREMAL ANNUAL REGIME OF THE SERIES
%%%%%%%%
if exist(['figures\PaperThreshold'],'dir')~=7,
    mkdir(['figures\PaperThreshold']) % Create a directory for figures if it doesn't exist
end

ejemplo='PRCP_ASN00021043';

% If it exists, read the file
if strcmp(ejemplo,'PRCP_ASN00021043'),
    % If the file exists, proceed to read it
    data=readtable([actualpath '\data\PRCP_ASN00021043.csv'],'TreatAsEmpty',{'NA'},'format','%s%f');
    % Before doing anything, remove null or NaN records
    % Verify the second column (precipitation)
    rowsToDelete = isnan(data{:, 2});

    % Remove rows that meet the condition
    data(rowsToDelete, :) = [];
    
    % Assign the data to the pluviometer structure
    pluviometros.nombre = 'PRCP_ASN00021043';
    pluviometros.datenum =datenum(data{:,1});
    pluviometros.fechas =data{:,1};
    pluviometros.data =data{:,2}/10; % Convert units if necessary
    
    % Calculate the minimum time difference between consecutive records
    % This is assumed to be in days
    min_time_diff = min(diff(pluviometros.datenum));
    
    % Define the minimum distance between peaks in temporal units
    % For example, the minimum distance will be twice the minimum distance between records
    min_peak_distance = 2 * min_time_diff;
    
    % Save the execution with 'n' in the file workspaceMC_POT_AM_n.mat
% % %     load(['PRCP_ASN00021043.mat'])
    threshold0 = 2.0;
end

% Define an initial threshold
threshold=0.0; 
n0=10; % Minimum number of peaks for valid statistics
siglevel=0.01; % Significance level for statistical tests
% Extract all valid thresholds in order and their corresponding
% mean exceedances and weights
[pks_unicos_valid,excedencias_mean_valid,excedencias_weight_valid, pks, locs, autocorrelations] = threshold_peak_extraction(pluviometros.data,threshold,n0,min_peak_distance);

% Check independence
for i=1:size(autocorrelations,1),
    if autocorrelations(i,3)<siglevel,
        fprintf('Lag %d: Significant autocorrelation (%.4f), p-val (%.4f)\n', ...
                autocorrelations(i,1), autocorrelations(i,2), autocorrelations(i,3));
    else
        fprintf('Lag %d: No Significant autocorrelation (%.4f), p-val (%.4f)\n', ...
                autocorrelations(i,1), autocorrelations(i,2), autocorrelations(i,3));
    end
end

% Optional: Plot the studentized residuals
if graficos,
    % Create the figure
    fonsiz = 18;
    scrsz = get(0, 'ScreenSize');
    figure('Position', [1 1 scrsz(3) scrsz(4)]);
    

    % Plot the unique independent peaks vs. mean exceedances using scatter plot
    scatter(pks_unicos_valid, excedencias_mean_valid, 'k', 'filled');  % Black filled markers

    % Customize axis labels with LaTeX
    xlabel('Threshold $u$ (mm/d)', 'FontSize', fonsiz, 'Interpreter', 'latex');
    ylabel('Mean exceedance $\bar{e}(u)$', 'FontSize', fonsiz, 'Interpreter', 'latex');

    % Add a title with LaTeX formatting
    title('Mean Residual Life Plot', 'FontSize', fonsiz+2, 'Interpreter', 'latex');

    % Enable the grid
    grid on;

    % Customize the axes and set LaTeX fonts
    set(gca, 'FontSize', fonsiz, 'TickLabelInterpreter', 'latex');

    % Keep the plot clean
    hold off;

    saveas(gcf, ['figures\PaperThreshold\MRLP' ejemplo], 'png');
    
    % Convert datenum to datetime for proper formatting in years
    years_ = datetime(pluviometros.datenum, 'ConvertFrom', 'datenum');
    
    % Find the first and last date in the 'years' array
    start_date = min(years_);
    end_date = max(years_);

    % Compute the difference in years
    years_diff = years(end_date - start_date);  % Approximate number of days in a year

    % Display the result
    disp(['The time span is approximately ' num2str(years_diff) ' years.']);

    % Create the figure
    fonsiz = 18;
    scrsz = get(0, 'ScreenSize');
    figure('Position', [1 1 scrsz(3) scrsz(4)]);

    % Plot the precipitation time series in gray
    
    plot(years_, pluviometros.data, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    hold on;

    % Plot the independent peaks with black triangles
    plot(years_(locs), pluviometros.data(locs), 'vk', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

    % Customize axis labels with LaTeX
    xlabel('Year', 'FontSize', fonsiz, 'Interpreter', 'latex');
    ylabel('Precipitation (mm)', 'FontSize', fonsiz, 'Interpreter', 'latex');

    % Add a legend
    legend({'Precipitation series', 'Independent peaks'}, 'Location', 'best', 'FontSize', fonsiz, 'Interpreter', 'latex');

    % Enable the grid
    grid on;

    % Format the plot with LaTeX fonts
    set(gca, 'FontSize', fonsiz, 'TickLabelInterpreter', 'latex');

    % Keep the plot clean
    hold off;

    saveas(gcf, ['figures\PaperThreshold\Series' ejemplo], 'png');

    
end

%% Method using studentized residuals and the smoothing spline
if graficos,
    [threshold_val_SR,beta,fobj,r] = threshold_studentized_residuals(pks_unicos_valid, excedencias_mean_valid, excedencias_weight_valid, siglevel, true,['figures\PaperThreshold\' ejemplo],false);
else
    [threshold_val_SR,beta,fobj,r] = threshold_studentized_residuals(pks_unicos_valid, excedencias_mean_valid, excedencias_weight_valid, siglevel);
end

%% Langousis method
if graficos,
    threshold_val_MSE = threshold_MSE(pks_unicos_valid, excedencias_mean_valid, excedencias_weight_valid, n0, [], true,['figures\PaperThreshold\' ejemplo]);
else
    threshold_val_MSE = threshold_MSE(pks_unicos_valid, excedencias_mean_valid, excedencias_weight_valid, n0);
end

%% Anderson-Darling method
threshold_val_AD = threshold_AD(pluviometros.data, siglevel);
threshold_val_AD = threshold_AD(pks, siglevel);

%% Cramer-Von Mises method
threshold_val_CVM = threshold_CVM(pluviometros.data, siglevel);
threshold_val_CVM = threshold_CVM(pks, siglevel);

[threshold_val_SR,threshold_val_MSE,threshold_val_AD,threshold_val_CVM]
