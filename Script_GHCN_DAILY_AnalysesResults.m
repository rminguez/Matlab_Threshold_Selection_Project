% Clear workspace and close all figures
clear all;
close all;
clc;

% Create 'figures' directory if it doesn't exist
if ~exist('figures', 'dir')
    mkdir('figures');
end

%%
% Read the threshold results from the CSV file
data = readtable('threshold_results.csv');

data.Properties.VariableNames

% Extract relevant columns for boxplot
% Combine all threshold columns into one matrix for easier plotting
threshold_data = [data.Langousis_Min1(data.Siglevel==0.05), data.Langousis_Min2(data.Siglevel==0.05), data.Langousis_Min3(data.Siglevel==0.05), ...
                  data.Studentized_Residuals(data.Siglevel==0.01),data.Studentized_Residuals(data.Siglevel==0.05),...
                  data.Anderson_Darling(data.Siglevel==0.01),data.Anderson_Darling(data.Siglevel==0.05),...
                  data.CVM(data.Siglevel==0.01),data.CVM(data.Siglevel==0.05)];

% Define labels for the methods
method_labels = {'Langousis Min 1', 'Langousis Min 2', 'Langousis Min 3', ...
                 'Studentized Residuals 1%', 'Studentized Residuals 5%', 'Anderson Darling 1%',...
                 'Anderson Darling 5%','Cramer-Von Mises 1%','Cramer-Von Mises 5%'};

% Create a boxplot with custom colors
fonsiz = 18;
scrsz = get(0, 'ScreenSize');
figure('Position', [1 1 scrsz(3) scrsz(4)]);
h = boxplot(threshold_data, 'Labels', method_labels, 'Colors', 'k', 'Widths', 0.5);

% Customize the boxplot appearance
set(h, 'LineWidth', 1.5);  % Thicker box edges

% Set gray color for the box face
colors = [0.8, 0.8, 0.8];  % Light gray color
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.5);  % Light gray fill with transparency
end

% Customize outliers to be black
outliers = findobj(gca, 'Tag', 'Outliers');
set(outliers, 'MarkerEdgeColor', 'k');  % Set outlier color to black

% Manually rotate x-axis labels by adjusting their properties
ax = gca;  % Get the current axis
xtick = get(ax, 'XTick');  % Get current x-tick positions
xticklabel = get(ax, 'XTickLabel');  % Get current x-tick labels
set(ax, 'XTickLabel', []);  % Remove original x-tick labels

% Adjust the y-position of the labels (move upwards slightly to make sure they fit in the figure)
y_offset = ax.YLim(1) - 0.05 * range(ax.YLim);  % Adjust the offset as needed

% Create rotated labels at the correct positions
text(xtick, repmat(y_offset, size(xtick)), xticklabel, ...
    'HorizontalAlignment', 'right', 'Rotation', 45, 'FontSize', 10, 'FontName', 'Times', 'Interpreter', 'latex');

% Set title, labels, and grid, using LaTeX fonts for all text
title('Threshold Selection Comparison Across Methods', 'FontWeight', 'Bold', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Selected Threshold (mm)', 'FontSize', 12, 'Interpreter', 'latex');
xlabel('', 'FontSize', 12, 'Interpreter', 'latex');
grid on;

% Customize the font and axis
set(gca, 'FontSize', 12, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');  % Use LaTeX for tick labels
set(gca, 'XColor', 'k', 'YColor', 'k');  % Black axis colors

% Adjust the figure to ensure the labels are visible
set(gca, 'Position', [0.13 0.25 0.75 0.7]);  % Adjust the axis position

% Save the figure in the 'figures' folder
saveas(gcf, 'figures/threshold_boxplot_gray_black_outliers_rotated.png');


%%
% Extract relevant columns for boxplot
% Combine all threshold columns into one matrix for easier plotting
time_data = [data.Langousis_Time(data.Siglevel==0.05), ...
                  data.Studentized_Residuals_Time(data.Siglevel==0.01),data.Studentized_Residuals_Time(data.Siglevel==0.05),...
                  data.Anderson_Darling_Time(data.Siglevel==0.01),data.Anderson_Darling_Time(data.Siglevel==0.05),...
                  data.CVM_Time(data.Siglevel==0.01),data.CVM_Time(data.Siglevel==0.05)];
              
% Define labels for the methods
time_labels = {'Langousis', ...
                 'Studentized Residuals 1%', 'Studentized Residuals 5%', 'Anderson Darling 1%',...
                 'Anderson Darling 5%','Cramer-Von Mises 1%','Cramer-Von Mises 5%'};

% Create a boxplot with custom colors
fonsiz = 18;
scrsz = get(0, 'ScreenSize');
figure('Position', [1 1 scrsz(3) scrsz(4)]);
h = boxplot(time_data, 'Labels', time_labels, 'Colors', 'k', 'Widths', 0.5);

% Customize the boxplot appearance
set(h, 'LineWidth', 1.5);  % Thicker box edges

% Set gray color for the box face
colors = [0.8, 0.8, 0.8];  % Light gray color
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.5);  % Light gray fill with transparency
end

% Customize outliers to be black
outliers = findobj(gca, 'Tag', 'Outliers');
set(outliers, 'MarkerEdgeColor', 'k');  % Set outlier color to black

% Manually rotate x-axis labels by adjusting their properties
ax = gca;  % Get the current axis
xtick = get(ax, 'XTick');  % Get current x-tick positions
xticklabel = get(ax, 'XTickLabel');  % Get current x-tick labels
set(ax, 'XTickLabel', []);  % Remove original x-tick labels

% Adjust the y-position of the labels (move upwards slightly to make sure they fit in the figure)
y_offset = ax.YLim(1) - 0.05 * range(ax.YLim);  % Adjust the offset as needed

% Create rotated labels at the correct positions
text(xtick, repmat(y_offset, size(xtick)), xticklabel, ...
    'HorizontalAlignment', 'right', 'Rotation', 45, 'FontSize', 10, 'FontName', 'Times', 'Interpreter', 'latex');

% Set title, labels, and grid, using LaTeX fonts for all text
title('Threshold Selection CPU  Time Comparison Across Methods', 'FontWeight', 'Bold', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Time (s.)', 'FontSize', 12, 'Interpreter', 'latex');
xlabel('', 'FontSize', 12, 'Interpreter', 'latex');
grid on;

% Customize the font and axis
set(gca, 'FontSize', 12, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');  % Use LaTeX for tick labels
set(gca, 'XColor', 'k', 'YColor', 'k');  % Black axis colors

% Adjust the figure to ensure the labels are visible
set(gca, 'Position', [0.13 0.25 0.75 0.7]);  % Adjust the axis position

% Save the figure in the 'figures' folder
saveas(gcf, 'figures/time_boxplot_gray_black_outliers_rotated.png');
              
%%

% Combine the threshold data into a matrix with each column representing a method
threshold_data = [data.Langousis_Min1(data.Siglevel==0.05), ...
                  data.Langousis_Min2(data.Siglevel==0.05), ...
                  data.Langousis_Min3(data.Siglevel==0.05), ...
                  data.Studentized_Residuals(data.Siglevel==0.01), ...
                  data.Studentized_Residuals(data.Siglevel==0.05), ...
                  data.Anderson_Darling(data.Siglevel==0.01), ...
                  data.Anderson_Darling(data.Siglevel==0.05), ...
                  data.CVM(data.Siglevel==0.01), ...
                  data.CVM(data.Siglevel==0.05)];

% Method labels for the table
method_labels = {'Langousis Min1', 'Langousis Min2', 'Langousis Min3', ...
                 'Studentized Residuals 1%', 'Studentized Residuals 5%', ...
                 'Anderson-Darling 1%', 'Anderson-Darling 5%', ...
                 'CVM 1%', 'CVM 5%'};

% Preallocate arrays to store summary statistics
mean_vals = zeros(length(method_labels), 1);
median_vals = zeros(length(method_labels), 1);
std_vals = zeros(length(method_labels), 1);
min_vals = zeros(length(method_labels), 1);
max_vals = zeros(length(method_labels), 1);
Q1_vals = zeros(length(method_labels), 1);
Q3_vals = zeros(length(method_labels), 1);
Critlo_vals = zeros(length(method_labels), 1);
Critup_vals = zeros(length(method_labels), 1);

% Calculate statistics for each method
for i = 1:length(method_labels)
    data_column = threshold_data(:, i);
    
    % Remove NaN values before calculating statistics
    data_clean = data_column(~isnan(data_column));
    
    % Compute the summary statistics
    mean_vals(i) = mean(data_clean);
    median_vals(i) = median(data_clean);
    std_vals(i) = std(data_clean);
    min_vals(i) = min(data_clean);
    max_vals(i) = max(data_clean);
    Q1_vals(i) = quantile(data_clean,0.25);
    Q3_vals(i) = quantile(data_clean,0.75);
    Critlo_vals(i) = max(0,Q1_vals(i) -1.5*(Q3_vals(i) -Q1_vals(i) ));
    Critup_vals(i) = Q3_vals(i) +1.5*(Q3_vals(i) -Q1_vals(i) );
end

% Create a table to store the results
summary_stats = table(method_labels', mean_vals, median_vals, std_vals, min_vals, max_vals,Q1_vals,Q3_vals,Critlo_vals,Critup_vals, ...
    'VariableNames', {'Method', 'Mean', 'Median', 'StdDev', 'Min', 'Max','Q1','Q3','Crit_lo','Crit_up'});

% Display the summary statistics
disp(summary_stats);

%%
% Extract the threshold data for the different methods
threshold_data = [data.Langousis_Min1(data.Siglevel==0.05); 
                  data.Langousis_Min2(data.Siglevel==0.05); 
                  data.Langousis_Min3(data.Siglevel==0.05); 
                  data.Studentized_Residuals(data.Siglevel==0.01); 
                  data.Studentized_Residuals(data.Siglevel==0.05); 
                  data.Anderson_Darling(data.Siglevel==0.01); 
                  data.Anderson_Darling(data.Siglevel==0.05); 
                  data.CVM(data.Siglevel==0.01); 
                  data.CVM(data.Siglevel==0.05)];

% Create a grouping variable for the method names
methods = [repmat({'Langousis_Min1'}, sum(data.Siglevel==0.05), 1);
           repmat({'Langousis_Min2'}, sum(data.Siglevel==0.05), 1);
           repmat({'Langousis_Min3'}, sum(data.Siglevel==0.05), 1);
           repmat({'Studentized_Residuals_1%'}, sum(data.Siglevel==0.01), 1);
           repmat({'Studentized_Residuals_5%'}, sum(data.Siglevel==0.05), 1);
           repmat({'Anderson_Darling_1%'}, sum(data.Siglevel==0.01), 1);
           repmat({'Anderson_Darling_5%'}, sum(data.Siglevel==0.05), 1);
           repmat({'CVM_1%'}, sum(data.Siglevel==0.01), 1);
           repmat({'CVM_5%'}, sum(data.Siglevel==0.05), 1)];

% Perform one-way ANOVA
[p, tbl, stats] = anova1(threshold_data, methods);

% Display the ANOVA table
disp('ANOVA Table:');
disp(tbl);

% Optional: Perform multiple comparison tests
multcompare(stats);

%%  Make some specific tests about the means and distributions
% Extract the relevant data for Langousis Min1 and Studentized Residuals 1%
langousis_min1 = data.Langousis_Min1(data.Siglevel == 0.01);
studentized_residuals_1 = data.Studentized_Residuals(data.Siglevel == 0.01);

% Remove NaN values from both datasets
langousis_min1_clean = langousis_min1(~isnan(langousis_min1));
studentized_residuals_1_clean = studentized_residuals_1(~isnan(studentized_residuals_1));

% Perform a two-sample t-test to compare the means
[H, P, CI, STATS] = ttest2(langousis_min1_clean, studentized_residuals_1_clean);

% Display the results
disp('T-test Results:');
disp(['Hypothesis Test Result (H): ', num2str(H)]);
disp(['P-value (P): ', num2str(P)]);
disp(['Confidence Interval (CI): ', num2str(CI')]);
disp(['Test Statistics (TSTAT): ', num2str(STATS.tstat)]);
disp(['Degrees of Freedom (DF): ', num2str(STATS.df)]);
disp(['Standard Deviation of the Difference (SD): ', num2str(STATS.sd)]);


% Perform the Kolmogorov-Smirnov test
[H, P, KSSTAT] = kstest2(langousis_min1_clean, studentized_residuals_1_clean);

% Display the results
disp('Kolmogorov-Smirnov Test Results:');
disp(['Hypothesis Test Result (H): ', num2str(H)]);
disp(['P-value (P): ', num2str(P)]);
disp(['KS Statistic: ', num2str(KSSTAT)]);


% Extract the relevant data for Langousis Min1 and Studentized Residuals 5%
langousis_min1 = data.Langousis_Min1(data.Siglevel == 0.05);
studentized_residuals_1 = data.Studentized_Residuals(data.Siglevel == 0.05);

% Remove NaN values from both datasets
langousis_min1_clean = langousis_min1(~isnan(langousis_min1));
studentized_residuals_1_clean = studentized_residuals_1(~isnan(studentized_residuals_1));

% Perform a two-sample t-test to compare the means
[H, P, CI, STATS] = ttest2(langousis_min1_clean, studentized_residuals_1_clean);

% Display the results
disp('T-test Results:');
disp(['Hypothesis Test Result (H): ', num2str(H)]);
disp(['P-value (P): ', num2str(P)]);
disp(['Confidence Interval (CI): ', num2str(CI')]);
disp(['Test Statistics (TSTAT): ', num2str(STATS.tstat)]);
disp(['Degrees of Freedom (DF): ', num2str(STATS.df)]);
disp(['Standard Deviation of the Difference (SD): ', num2str(STATS.sd)]);


% Perform the Kolmogorov-Smirnov test
[H, P, KSSTAT] = kstest2(langousis_min1_clean, studentized_residuals_1_clean);

% Display the results
disp('Kolmogorov-Smirnov Test Results:');
disp(['Hypothesis Test Result (H): ', num2str(H)]);
disp(['P-value (P): ', num2str(P)]);
disp(['KS Statistic: ', num2str(KSSTAT)]);
