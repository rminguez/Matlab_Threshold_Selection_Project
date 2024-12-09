% Clear workspace and close all figures
clear all;
close all;
clc;

% Define data and figures directories
data_folder = 'data';
figures_folder = 'figures';

% Ensure directories exist
if ~exist(data_folder, 'dir')
    mkdir(data_folder);
end
if ~exist(figures_folder, 'dir')
    mkdir(figures_folder);
end

%% Load simulation data
% Read the simulation results from the CSV file
data = readtable('data\simulation_results_with_rounding_and_runs.csv');

% Display variable names to check structure
disp(data.Properties.VariableNames);

%% Boxplot of selected thresholds by method and significance level
% Extract relevant threshold data for boxplot
threshold_data = [data.Langousis_Min1(data.SigLevel == 0.01), data.Langousis_Min2(data.SigLevel == 0.01), data.Langousis_Min3(data.SigLevel == 0.01), ...
                  data.Studentized_Residuals(data.SigLevel == 0.01), ...
                  data.Studentized_Residuals(data.SigLevel == 0.05), ...
                  data.Anderson_Darling(data.SigLevel == 0.01), ...
                  data.Anderson_Darling(data.SigLevel == 0.05), ...
                  data.Cramer_Von_Mises(data.SigLevel == 0.01), ...
                  data.Cramer_Von_Mises(data.SigLevel == 0.05)] - repmat(data.Real_Threshold(data.SigLevel == 0.01), 1, 9);

% Define method labels with significance levels
method_labels = {
                 'Langousis_Min1', 'Langousis_Min2', 'Langousis_Min3', ...
                 'Studentized_Residuals_1','Studentized_Residuals_5', ...
                 'Anderson_Darling_1','Anderson_Darling_5', ...
                 'Cramer_Von_Mises_1','Cramer_Von_Mises_5'};


% Define method labels
plot_labels = {'L1', 'L2', 'L3', ...
                 'SR 1\%', 'SR 5\%', 'AD 1\%', ...
                 'AD 5\%', 'CVM 1\%', 'CVM 5\%'};

% Compute percentages of NaN values and absolute error < 1
num_methods = size(threshold_data, 2);
num_cases = size(threshold_data, 1);

nan_percentage = zeros(1, num_methods);
error_within_one_percentage = zeros(1, num_methods);

for i = 1:num_methods
    % Calculate percentage of NaN values
    nan_percentage(i) = sum(isnan(threshold_data(:, i))) / num_cases * 100;
    
    % Calculate percentage of cases with |error| < 1
    valid_data = threshold_data(~isnan(threshold_data(:, i)), i); % Exclude NaN values
    error_within_one_percentage(i) = sum(abs(valid_data) < 1) / length(valid_data) * 100;
end

% Display the percentages
disp('Percentage of NaN values by method:');
disp(array2table(nan_percentage, 'VariableNames', method_labels));

disp('Percentage of |error| < 1 by method:');
disp(array2table(error_within_one_percentage, 'VariableNames', method_labels));

% Create a boxplot for thresholds
fonsiz = 24;
scrsz = get(0, 'ScreenSize');
figure('Position', [1 1 scrsz(3) scrsz(4)]);
h = boxplot(threshold_data, 'Labels', plot_labels, 'Colors', 'k', 'Widths', 0.5);

% Customize boxplot appearance
set(h, 'LineWidth', 1.5);
colors = [0.8, 0.8, 0.8]; % Light gray
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.5);
end

% Customize outliers to be black
outliers = findobj(gca, 'Tag', 'Outliers');
set(outliers, 'MarkerEdgeColor', 'k');  % Set outlier color to black

% Customize labels and grid
title('Threshold Selection Comparison Across Methods and Significance Levels', ...
      'FontWeight', 'Bold', 'FontSize', fonsiz + 2, 'Interpreter', 'latex');
ylabel('Selected Threshold', 'FontSize', fonsiz, 'Interpreter', 'latex');
xlabel('Method and Significance Level', 'FontSize', fonsiz, 'Interpreter', 'latex');
grid on;

% Customize tick labels
set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

% Add percentage text on top of each boxplot
xticks = 1:length(plot_labels); % X-axis positions for each method
for i = 1:length(xticks)
    percentage_text = sprintf('NaN: %.1f%%\n<1: %.1f%%', nan_percentage(i), error_within_one_percentage(i));
    
    text(xticks(i), max(threshold_data(:)) - 0.5, percentage_text, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', fonsiz - 6, 'Interpreter', 'latex', 'FontName', 'Times');
end

% Save the figure
saveas(gcf, [figures_folder '\threshold_boxplot_with_percentages.png']);



%% Boxplot of selected thresholds by method and significance level for different rounding values
% Get unique rounding values
rounding_values = unique(data.Rounding);

method_labels = {'Studentized Residuals', ...
                 'Langousis Min1', 'Langousis Min2', 'Langousis Min3', ...
                 'Anderson-Darling', 'Cramer-Von Mises'};

% Define method labels
plot_labels = {'L1', 'L2', 'L3', ...
                 'SR 1\%', 'SR 5\%', 'AD 1\%', ...
                 'AD 5\%', 'CVM 1\%', 'CVM 5\%'};

% Font size for plots
fonsiz = 24; 
scrsz = get(0, 'ScreenSize'); % Get screen size for figure positioning

for r = 1:length(rounding_values)
    % Filter data for the current rounding value
    current_rounding = rounding_values(r);
    filtered_data = data(data.Rounding == current_rounding, :);
    
    % Compute deviations from the real threshold
    threshold_data = [filtered_data.Langousis_Min1(filtered_data.SigLevel == 0.01), filtered_data.Langousis_Min2(filtered_data.SigLevel == 0.01), filtered_data.Langousis_Min3(filtered_data.SigLevel == 0.01), ...
                      filtered_data.Studentized_Residuals(filtered_data.SigLevel == 0.01), ...
                      filtered_data.Studentized_Residuals(filtered_data.SigLevel == 0.05), ...
                      filtered_data.Anderson_Darling(filtered_data.SigLevel == 0.01), ...
                      filtered_data.Anderson_Darling(filtered_data.SigLevel == 0.05), ...
                      filtered_data.Cramer_Von_Mises(filtered_data.SigLevel == 0.01), ...
                      filtered_data.Cramer_Von_Mises(filtered_data.SigLevel == 0.05)] - ...
                      repmat(filtered_data.Real_Threshold(filtered_data.SigLevel == 0.01), 1, 9);

    % Compute NaN percentages and error percentages
    nan_percentages = sum(isnan(threshold_data)) / size(threshold_data, 1) * 100;
    error_within_one_percentage = sum(abs(threshold_data) < 1) / size(threshold_data, 1) * 100;

    % Create a boxplot for thresholds
    figure('Position', [1 1 scrsz(3) scrsz(4)]);
    h = boxplot(threshold_data, 'Labels', plot_labels, 'Colors', 'k', 'Widths', 0.5);

    % Customize boxplot appearance
    set(h, 'LineWidth', 1.5);
    colors = [0.8, 0.8, 0.8]; % Light gray
    boxes = findobj(gca, 'Tag', 'Box');
    for j = 1:length(boxes)
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.5);
    end

    % Add percentage text on top of each boxplot
    xticks = 1:length(plot_labels); % X-axis positions for each method
    for i = 1:length(xticks)
        percentage_text = sprintf('NaN: %.1f%%\n<1: %.1f%%', nan_percentages(i), error_within_one_percentage(i));
        text(xticks(i), max(threshold_data(:)) - 0.5, percentage_text, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', fonsiz - 6, 'Interpreter', 'latex', 'FontName', 'Times');
    end

    % Customize labels and grid
    title(sprintf('Threshold Selection Comparison (Rounding = %d)', current_rounding), ...
          'FontWeight', 'Bold', 'FontSize', fonsiz + 2, 'Interpreter', 'latex');
    ylabel('Deviation from Real Threshold', 'FontSize', fonsiz, 'Interpreter', 'latex');
    xlabel('Method and Significance Level', 'FontSize', fonsiz, 'Interpreter', 'latex');
    grid on;
    set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

    % Save the figure
    saveas(gcf, [figures_folder sprintf('\\threshold_boxplot_rounding_%d.png', current_rounding)]);

end


%% Boxplot of selected thresholds by method and significance level for different k values
% Boxplot of selected thresholds by method and significance level for different k values
% Get unique k values
k_values = unique(data.k);

% Define plot labels
plot_labels = {'L1', 'L2', 'L3', ...
               'SR 1\%', 'SR 5\%', 'AD 1\%', ...
               'AD 5\%', 'CVM 1\%', 'CVM 5\%'};

% Font size for plots
fonsiz = 24; 

% Get screen size for figure positioning
scrsz = get(0, 'ScreenSize');

% Loop through each k value and create a boxplot
for k_idx = 1:length(k_values)
    % Filter data for the current k value
    current_k = k_values(k_idx);
    filtered_data = data(data.k == current_k, :);
    
    % Compute deviations from the real threshold
    threshold_data = [filtered_data.Langousis_Min1(filtered_data.SigLevel == 0.01), filtered_data.Langousis_Min2(filtered_data.SigLevel == 0.01), filtered_data.Langousis_Min3(filtered_data.SigLevel == 0.01), ...
                      filtered_data.Studentized_Residuals(filtered_data.SigLevel == 0.01), ...
                      filtered_data.Studentized_Residuals(filtered_data.SigLevel == 0.05), ...
                      filtered_data.Anderson_Darling(filtered_data.SigLevel == 0.01), ...
                      filtered_data.Anderson_Darling(filtered_data.SigLevel == 0.05), ...
                      filtered_data.Cramer_Von_Mises(filtered_data.SigLevel == 0.01), ...
                      filtered_data.Cramer_Von_Mises(filtered_data.SigLevel == 0.05)] - ...
                      repmat(filtered_data.Real_Threshold(filtered_data.SigLevel == 0.01), 1, 9);
    
    % Calculate NaN percentages and error percentages
    nan_percentages = sum(isnan(threshold_data)) / size(threshold_data, 1) * 100;
    error_within_one_percentage = sum(abs(threshold_data) < 1) / size(threshold_data, 1) * 100;

    % Create a boxplot for thresholds
    figure('Position', [1 1 scrsz(3) scrsz(4)]);
    h = boxplot(threshold_data, 'Labels', plot_labels, 'Colors', 'k', 'Widths', 0.5);

    % Customize boxplot appearance
    set(h, 'LineWidth', 1.5);
    colors = [0.8, 0.8, 0.8]; % Light gray
    boxes = findobj(gca, 'Tag', 'Box');
    for j = 1:length(boxes)
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.5);
    end

    % Add percentage text for NaN and absolute error < 1
    xticks = 1:length(plot_labels); % X-axis positions for each method
    for i = 1:length(xticks)
        percentage_text = sprintf('NaN: %.1f%%\n<1: %.1f%%', nan_percentages(i), error_within_one_percentage(i));
        text(xticks(i), max(threshold_data(:)) - 0.5, percentage_text, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', fonsiz - 6, 'Interpreter', 'latex', 'FontName', 'Times');
    end

    % Customize labels and grid
    title(sprintf('Threshold Selection Comparison ($\\xi = %.2f$)', current_k), ...
          'FontWeight', 'Bold', 'FontSize', fonsiz + 2, 'Interpreter', 'latex');
    ylabel('Deviation from Real Threshold', 'FontSize', fonsiz, 'Interpreter', 'latex');
    xlabel('Method and Significance Level', 'FontSize', fonsiz, 'Interpreter', 'latex');
    grid on;
    set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

    % Save the figure
    saveas(gcf, [figures_folder sprintf('\\threshold_boxplot_k_%.2f.png', current_k)]);

end



%% Boxplot of computational times
% Extract computational time data for boxplot
time_data = [data.SR_Time, data.Langousis_Time, ...
             data.AD_Time, data.CVM_Time];

% Define method labels
time_labels = {'Studentized Residuals', 'Langousis', ...
               'Anderson-Darling', 'Cramer-Von Mises'};

% Font size for plots
fonsiz = 24; 
           
% Create a boxplot for computational times
figure('Position', [1 1 scrsz(3) scrsz(4)]);
h = boxplot(time_data, 'Labels', time_labels, 'Colors', 'k', 'Widths', 0.5);

% Customize boxplot appearance
set(h, 'LineWidth', 1.5);
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.5);
end

% Customize labels and grid
title('Computational Time Comparison Across Methods', 'FontWeight', 'Bold', 'FontSize', fonsiz+2, 'Interpreter', 'latex');
ylabel('Time (s)', 'FontSize', fonsiz, 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

% Save the figure
saveas(gcf, [figures_folder '\\time_boxplot_with_runs.png']);

%% Summary statistics for thresholds
methods = {'Studentized Residuals', 'Langousis Min1', 'Langousis Min2', ...
           'Langousis Min3', 'Anderson-Darling', 'Cramer-Von Mises'};

mean_vals = zeros(length(methods), 1);
median_vals = zeros(length(methods), 1);
std_vals = zeros(length(methods), 1);

for i = 1:length(methods)
    % Extract method-specific data
    method_data = threshold_data(:, i);
    method_data_clean = method_data(~isnan(method_data)); % Remove NaNs

    % Compute statistics
    mean_vals(i) = mean(method_data_clean);
    median_vals(i) = median(method_data_clean);
    std_vals(i) = std(method_data_clean);
end

% Create a table with summary statistics
summary_stats = table(methods', mean_vals, median_vals, std_vals, ...
    'VariableNames', {'Method', 'Mean', 'Median', 'StdDev'});

% Display summary statistics
disp('Summary Statistics for Thresholds:');
disp(summary_stats);

%% Statistical tests: Studentized Residuals vs Langousis Min1
% Extract data for the comparison
studentized_residuals = data.Studentized_Residuals;
langousis_min1 = data.Langousis_Min1;

% Remove NaN values
studentized_residuals_clean = studentized_residuals(~isnan(studentized_residuals));
langousis_min1_clean = langousis_min1(~isnan(langousis_min1));

% Perform two-sample t-test
[H, P, CI, STATS] = ttest2(studentized_residuals_clean, langousis_min1_clean);

% Display t-test results
disp('T-test Results: Studentized Residuals vs Langousis Min1');
disp(['Hypothesis Test Result (H): ', num2str(H)]);
disp(['P-value (P): ', num2str(P)]);
disp(['Confidence Interval (CI): ', num2str(CI')]);
disp(['Test Statistics (TSTAT): ', num2str(STATS.tstat)]);

%% Save summary statistics to CSV
writetable(summary_stats, 'data\summary_statistics_with_runs.csv');
disp('Summary statistics saved to summary_statistics_with_runs.csv');

%% ANOVA ANALYSIS 
% 
% Compute the difference between the estimated thresholds and the real threshold
methods = {'Studentized_Residuals', 'Langousis_Min1', 'Langousis_Min2', 'Langousis_Min3', ...
           'Anderson_Darling', 'Cramer_Von_Mises'};
threshold_diff = [];

% Reshape data to long format for ANOVA
for i = 1:length(methods)
    method_data = data.(methods{i}) - data.Real_Threshold;
    threshold_diff = [threshold_diff; method_data];
end

% Create grouping variables
method_group = repelem(methods, height(data))';
rounding_group = repelem(data.Rounding, length(methods));

% Prepare table for ANOVA
anova_data = table(threshold_diff, method_group, rounding_group, ...
    'VariableNames', {'ThresholdDiff', 'Method', 'Rounding'});

% Perform two-way ANOVA
[p, tbl, stats] = anovan(anova_data.ThresholdDiff, ...
    {anova_data.Method, anova_data.Rounding}, ...
    'model', 'interaction', 'varnames', {'Method', 'Rounding'});

% Display the ANOVA table
disp('ANOVA Table:');
disp(tbl);

% Multiple comparisons for Method
figure;
[c_method, ~, ~, ~] = multcompare(stats, 'Dimension', 1);
title('Multiple Comparisons: Method', 'Interpreter', 'latex');
ylabel('Threshold Difference', 'FontSize', 12, 'Interpreter', 'latex');

% Display results
disp('Pairwise comparison results:');
disp(c_method);

% Multiple comparisons for Rounding
figure;
[c_rounding, ~, ~, ~] = multcompare(stats, 'Dimension', 2);
title('Multiple Comparisons: Rounding', 'Interpreter', 'latex');
ylabel('Threshold Difference', 'FontSize', 12, 'Interpreter', 'latex');


%%
% Prepare data for grouped boxplot
methods = {'Studentized_Residuals', 'Langousis_Min1', 'Langousis_Min2', 'Langousis_Min3', ...
           'Anderson_Darling', 'Cramer_Von_Mises'};

lmethods = {'SR', 'L1', 'L2', 'L3', 'AD', 'CV'};       
threshold_diff = [];
method_group = [];
rounding_group = [];

% Reshape data
for i = 1:length(methods)
    method_data = filtered_data.(methods{i}) - filtered_data.Real_Threshold;
    threshold_diff = [threshold_diff; method_data];
    method_group = [method_group; repelem({lmethods{i}}, length(method_data))'];
    rounding_group = [rounding_group; filtered_data.Rounding];
end

% Convert grouping variables to categorical
method_group = categorical(method_group);
rounding_group = categorical(rounding_group);

% Create the grouped boxplot
fonsiz = 20;
scrsz = get(0, 'ScreenSize');
figure('Position', [1 1 scrsz(3) scrsz(4)]);
boxplot(threshold_diff, {method_group, rounding_group}, ...
    'Colors', 'k', 'FactorSeparator', 1, 'Widths', 0.5);

% Customize boxplot appearance
set(findobj(gca, 'type', 'line'), 'linew', 1.5); % Thicker lines
xlabel('Method and Rounding', 'FontSize', fonsiz, 'Interpreter', 'latex');
ylabel('Threshold Difference', 'FontSize', fonsiz, 'Interpreter', 'latex');
title('Threshold Difference by Method and Rounding for k=0.15', 'FontWeight', 'Bold', ...
    'FontSize', fonsiz, 'Interpreter', 'latex');
grid on;

% Adjust x-axis tick labels manually
ax = gca;
xticks = get(ax, 'XTick');
xtick_labels = get(ax, 'XTickLabel');
set(ax, 'XTickLabel', []); % Remove default labels
xtick_labels = cellstr(xtick_labels); % Convert labels to cell array
y_offset = ax.YLim(1) - 0.05 * range(ax.YLim); % Adjust position for labels

% Add rotated labels
for i = 1:length(xticks)
    text(xticks(i), y_offset, xtick_labels{i}, 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'Rotation', 45, 'FontSize', fonsiz, 'Interpreter', 'latex');
end

% Adjust spacing to avoid overlap
ax.XTickLabelRotation = 0; % Ensure rotation doesn't interfere
ax.XTick = 1:5:length(xticks); % Adjust tick intervals for readability
set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

% Save the figure
saveas(gcf, [figures_folder '\\threshold_diff_boxplot_by_method_and_rounding_fixed.png']);

 