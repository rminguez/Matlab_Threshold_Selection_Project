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
data = readtable('data\threshold_results.csv');

data.Properties.VariableNames

% Extract relevant columns for boxplot
% Combine all threshold columns into one matrix for easier plotting
threshold_data = [data.Langousis_Min1(data.Siglevel==0.05), data.Langousis_Min2(data.Siglevel==0.05), data.Langousis_Min3(data.Siglevel==0.05), ...
                  data.Studentized_Residuals(data.Siglevel==0.01),data.Studentized_Residuals(data.Siglevel==0.05),...
                  data.Anderson_Darling(data.Siglevel==0.01),data.Anderson_Darling(data.Siglevel==0.05),...
                  data.CVM(data.Siglevel==0.01),data.CVM(data.Siglevel==0.05)];

% Define labels for the methods
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
    error_within_one_percentage(i) = sum((valid_data>=2 & valid_data<=12) ) / length(valid_data) * 100;
end

% Display the percentages
disp('Percentage of NaN values by method:');
disp(array2table(nan_percentage, 'VariableNames', method_labels));

disp('Percentage of |error| < 1 by method:');
disp(array2table(error_within_one_percentage, 'VariableNames', method_labels));

% Create a boxplot with custom colors
fonsiz = 24;
scrsz = get(0, 'ScreenSize');
figure('Position', [1 1 scrsz(3) scrsz(4)]);
h = boxplot(threshold_data, 'Labels', plot_labels, 'Colors', 'k', 'Widths', 0.5);

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

% Set title, labels, and grid, using LaTeX fonts for all text
title('Threshold Selection Comparison Across Methods', 'FontWeight', 'Bold', 'FontSize', fonsiz+2, 'Interpreter', 'latex');
ylabel('Selected Threshold (mm)', 'FontSize', fonsiz, 'Interpreter', 'latex');
xlabel('Method','FontSize', fonsiz, 'Interpreter', 'latex');
grid on;

% Customize tick labels
set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

% Add percentage text on top of each boxplot
xticks = 1:length(plot_labels); % X-axis positions for each method
for i = 1:length(xticks)
    percentage_text = sprintf('NaN: %.1f%%\n%.1f%%', nan_percentage(i), error_within_one_percentage(i));
    
    text(xticks(i), max(threshold_data(:)) - 10, percentage_text, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', fonsiz - 6, 'Interpreter', 'latex', 'FontName', 'Times');
end

% Optional: Save the figure
saveas(gcf, 'figures\threshold_boxplot_gray_black_outliers_rotated.png');

%%
% Extract relevant columns for boxplot
% Combine all threshold columns into one matrix for easier plotting
time_data = [data.Langousis_Time(data.Siglevel==0.05), ...
                  data.Studentized_Residuals_Time(data.Siglevel==0.01),data.Studentized_Residuals_Time(data.Siglevel==0.05),...
                  data.Anderson_Darling_Time(data.Siglevel==0.01),data.Anderson_Darling_Time(data.Siglevel==0.05),...
                  data.CVM_Time(data.Siglevel==0.01),data.CVM_Time(data.Siglevel==0.05)];
              
% Define labels for the times
time_labels = {
                 'Langousis',  ...
                 'Studentized_Residuals_1','Studentized_Residuals_5', ...
                 'Anderson_Darling_1','Anderson_Darling_5', ...
                 'Cramer_Von_Mises_1','Cramer_Von_Mises_5'};             
% Define method labels
time_plot_labels = {'L', ...
                 'SR 1\%', 'SR 5\%', 'AD 1\%', ...
                 'AD 5\%', 'CVM 1\%', 'CVM 5\%'};
% Create a boxplot with custom colors
fonsiz = 24;
scrsz = get(0, 'ScreenSize');
figure('Position', [1 1 scrsz(3) scrsz(4)]);
h = boxplot(time_data, 'Labels', time_plot_labels, 'Colors', 'k', 'Widths', 0.5);

% Customize the boxplot appearance
set(h, 'LineWidth', 1.5);  % Thicker box edges

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

% Set title, labels, and grid, using LaTeX fonts for all text
title('Threshold Selection CPU  Time Comparison Across Methods', 'FontWeight', 'Bold', 'FontSize', fonsiz+2, 'Interpreter', 'latex');
ylabel('Time (s.)', 'FontSize', fonsiz, 'Interpreter', 'latex');
xlabel('', 'FontSize', fonsiz, 'Interpreter', 'latex');
grid on;

% Customize the font and axis
set(gca, 'FontSize', fonsiz,  'TickLabelInterpreter', 'latex');  % Use LaTeX for tick labels
set(gca, 'XColor', 'k', 'YColor', 'k');  % Black axis colors

% Optional: Save the figure
saveas(gcf, 'figures\time_boxplot_gray_black_outliers_rotated.png');
              
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
method_labels = {'L1', 'L2', 'L3', ...
                 'SR 1\%', 'SR 5\%', 'AD 1\%',...
                 'AD 5\%','CVM 1\%','CVM 5\%'};
             
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
percentage_vals = zeros(length(method_labels), 1);
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
    percentage_vals(i) = sum((data_clean>=2 & data_clean<=12) ) / length(data_clean) * 100;
end

% Create a table to store the results
summary_stats = table(method_labels', mean_vals, median_vals, std_vals, min_vals, max_vals,Q1_vals,Q3_vals,Critlo_vals,Critup_vals,percentage_vals, ...
    'VariableNames', {'Method', 'Mean', 'Median', 'StdDev', 'Min', 'Max','Q1','Q3','Crit_lo','Crit_up','Percentage'});

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


%% Plot the selected stations on a world map
% Extract latitude and longitude of selected stations
load('data\selected_stations_data.mat');

% Plot the world map
fonsiz = 24;
scrsz = get(0, 'ScreenSize');
figure('Position', [1 1 scrsz(3) scrsz(4)]);
ax_ = newplot;
worldmap('World'); % Create a world map
plotm(coastlat, coastlon, 'k'); % Plot coastlines

% Add the station locations
hold on;
scatterm(latitudes, longitudes, 50, 'k', '^', 'filled'); % Use red filled triangles for markers
title('Selected Precipitation Stations', 'FontSize', fonsiz, 'Interpreter', 'latex');

% Optional: Customize the map
setm(gca, 'FontSize', 12, 'FontName', 'Helvetica');
gridm('on'); % Add grid lines
mlabel('on'); % Add meridian labels
plabel('on'); % Add parallel labels

% Save the map
saveas(gcf, 'figures\WorldMap_SelectedStations.png');




%%  Scatterm for thresholds
%Load the station data
load('data\selected_stations_data.mat'); % Contains latitudes, longitudes


threshold_data = [data.Langousis_Min1(data.Siglevel==0.05)  data.Langousis_Min2(data.Siglevel==0.05)  data.Langousis_Min3(data.Siglevel==0.05)  ...
                  data.Studentized_Residuals(data.Siglevel==0.01) data.Studentized_Residuals(data.Siglevel==0.05) ...
                  data.Anderson_Darling(data.Siglevel==0.01) data.Anderson_Darling(data.Siglevel==0.05) ...
                  data.CVM(data.Siglevel==0.01) data.CVM(data.Siglevel==0.05)];

% Define method labels
method_labels = {'L1', 'L2', 'L3', ...
                 'SR 1\%', 'SR 5\%', ...
                 'AD 1\%', 'AD 5\%', ...
                 'CVM 1\%', 'CVM 5\%'};

latitudes = latitudes(idx);
longitudes = longitudes(idx);



% Ensure latitudes, longitudes, and thresholds are the same size
if size(latitudes, 1) ~= size(threshold_data, 1)
    error('Latitudes, longitudes, and threshold data must have the same number of rows.');
end

% Loop through each method to create a map
fonsiz = 24;
scrsz = get(0, 'ScreenSize');



for method_idx = 1:length(method_labels)
    % Extract threshold data for the current method
    thresholds = threshold_data(:, method_idx);
    
    % Remove entries with NaN thresholds
    valid_idx = ~isnan(thresholds);
    valid_latitudes = latitudes(valid_idx);
    valid_longitudes = longitudes(valid_idx);
    valid_thresholds = thresholds(valid_idx);
    
    % Plot the world map
    figure('Position', [1 1 scrsz(3) scrsz(4)]);
    worldmap('World'); % Create a world map
    load coastlines; % Load coastlines for reference
    plotm(coastlat, coastlon, 'k'); % Plot coastlines
    
    hold on;
     % Add a colorbar
    colormap('gray'); % Change colormap to jet for heatmap
    c = colorbar;
    c.Label.String = 'Threshold (mm)';
    c.Label.FontSize = fonsiz;
    
    % Add station locations with heat map coloring

    scatterm(valid_latitudes, valid_longitudes, 50, valid_thresholds, 'filled'); % Color-coded thresholds
    
   
    
    % Set title
    title(['Thresholds: ', method_labels{method_idx}], 'FontSize', fonsiz, 'Interpreter', 'latex');
    
    % Customize the map
    setm(gca, 'FontSize', 12, 'FontName', 'Helvetica');
    gridm('on'); % Add grid lines
    mlabel('on'); % Add meridian labels
    plabel('on'); % Add parallel labels
    
    % Save the map
    saveas(gcf, ['figures\WorldMap_Thresholds_' strrep(strrep(method_labels{method_idx}, ' ', '_'),'\%','') '.png']);
end

close all

%% Define the range of the color bar
color_range = [-50, 150];

% Total number of colors in the colormap
n_colors = 256;

% Position of 0 in the color scale
midpoint = abs(color_range(1)) / (abs(color_range(1)) + color_range(2));

% Number of colors for the negative and positive parts
n_negative = round(midpoint * n_colors);
n_positive = n_colors - n_negative;

% Generate colors for the negative range (from white to blue)
red_negative = linspace(1, 0, n_negative)'; % Red component decreases from 1 to 0
green_negative = linspace(1, 0, n_negative)'; % Green component decreases from 1 to 0
blue_negative = ones(n_negative, 1); % Blue component remains constant at 1
negative_colors = [red_negative, green_negative, blue_negative];

% Generate colors for the positive range (from white to red)
red_positive = ones(n_positive, 1); % Red component remains constant at 1
green_positive = linspace(1, 0, n_positive)'; % Green component decreases from 1 to 0
blue_positive = linspace(1, 0, n_positive)'; % Blue component decreases from 1 to 0
positive_colors = [red_positive, green_positive, blue_positive];

% Combine the negative and positive colors
custom_colormap = [flipud(negative_colors); positive_colors];

%% Contour plots
%Load the station data
load('data\selected_stations_data.mat'); % Contains latitudes, longitudes

% Threshold data and labels
threshold_data = [data.Langousis_Min1(data.Siglevel==0.05), ...
                  data.Langousis_Min2(data.Siglevel==0.05), ...
                  data.Studentized_Residuals(data.Siglevel==0.01), ...
                  data.Studentized_Residuals(data.Siglevel==0.05), ...
                  data.Anderson_Darling(data.Siglevel==0.01), ...
                  data.Anderson_Darling(data.Siglevel==0.05), ...
                  data.CVM(data.Siglevel==0.01), ...
                  data.CVM(data.Siglevel==0.05)];
method_labels = {'L1', 'L2', ...
                 'SR 1\%', 'SR 5\%', ...
                 'AD 1\%', 'AD 5\%', ...
                 'CVM 1\%', 'CVM 5\%'};

% Valid locations
valid_latitudes = latitudes(idx);
valid_longitudes = longitudes(idx);





% Loop through each method to create a contour map
fonsiz = 24;
scrsz = get(0, 'ScreenSize');
grid_spacing = 1; % Grid spacing in degrees
max_distance = 2; % Maximum distance in degrees for interpolation

for method_idx = 1:length(method_labels)
    % Extract thresholds for the current method
    valid_thresholds = threshold_data(:, method_idx);
    
    
   
    
    % Create a grid of latitudes and longitudes
    [grid_lon, grid_lat] = meshgrid(-180:grid_spacing:180, -90:grid_spacing:90);
    
    % Interpolate threshold data
    F = scatteredInterpolant(valid_longitudes, valid_latitudes, valid_thresholds, 'linear', 'none');
    grid_values = F(grid_lon, grid_lat);
    
    % Calculate distances to nearest stations
    station_coords = [valid_longitudes, valid_latitudes];
    grid_coords = [grid_lon(:), grid_lat(:)];
    distances = pdist2(station_coords, grid_coords);
    min_distances = min(distances, [], 1);
    clear distances
    min_distances = reshape(min_distances, size(grid_lon));
    
    % Mask out regions where distance exceeds max_distance
    grid_values(min_distances > max_distance) = NaN;
    
    % Plot the world map
    figure('Position', [1 1 scrsz(3) scrsz(4)]);
    worldmap('World'); % Create a world map
    load coastlines; % Load coastlines for reference
    plotm(coastlat, coastlon, 'k'); % Plot coastlines
    
    % Add contour plot
    hold on;
%        scatterm(valid_latitudes, valid_longitudes, 50, valid_thresholds, 'filled'); % Color-coded thresholds
%   
    
    contourm(grid_lat, grid_lon, grid_values, 'LineColor', 'k'); % Contour lines
    pcolorm(grid_lat, grid_lon, grid_values); % Filled contour
    colormap(custom_colormap); % Jet colormap
    caxis([min(valid_thresholds), max(valid_thresholds)]); % Adjust color scale
    c = colorbar;
    c.Label.String = 'Threshold (mm)';
    c.Label.FontSize = fonsiz;
    
    % Set title
    title(['Interpolated Thresholds: ', method_labels{method_idx}], 'FontSize', fonsiz, 'Interpreter', 'latex');
    
    % Customize the map
    setm(gca, 'FontSize', 12, 'FontName', 'Helvetica');
    gridm('on'); % Add grid lines
    mlabel('on'); % Add meridian labels
    plabel('on'); % Add parallel labels
    
     
    set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex' );
   
    % Save the map
    saveas(gcf, ['figures\ContourMap_Thresholds_' strrep(strrep(method_labels{method_idx}, ' ', '_'), '\%', '') '.png']);
end
close all



%% Differences using scatterm and boxplots
% Load the station data
load('data\selected_stations_data.mat'); % Contains latitudes, longitudes

% Define threshold data columns for comparison
threshold_data = [data.Langousis_Min1(data.Siglevel==0.05), ...
                  data.Langousis_Min2(data.Siglevel==0.05), ...
                  data.Studentized_Residuals(data.Siglevel==0.01), ...
                  data.Studentized_Residuals(data.Siglevel==0.05), ...
                  data.Anderson_Darling(data.Siglevel==0.01), ...
                  data.Anderson_Darling(data.Siglevel==0.05), ...
                  data.CVM(data.Siglevel==0.01), ...
                  data.CVM(data.Siglevel==0.05)];

% Define method labels
method_labels = {'L2 - L1', ...
                 'SR 1\% - L1', ...
                 'SR 5\% - L1', ...
                 'AD 1\% - L1', ...
                 'AD 5\% - L1', ...
                 'CVM 1\% - L1', ...
                 'CVM 5\% - L1'};

% Extract Langousis Min1 as the baseline
langousis_min1 = threshold_data(:, 1);

% Calculate differences for each method
difference_data = [threshold_data(:, 2) - langousis_min1, ...
                   threshold_data(:, 3) - langousis_min1, ...
                   threshold_data(:, 4) - langousis_min1, ...
                   threshold_data(:, 5) - langousis_min1, ...
                   threshold_data(:, 6) - langousis_min1, ...
                   threshold_data(:, 7) - langousis_min1, ...
                   threshold_data(:, 8) - langousis_min1];

% Filter latitudes and longitudes
latitudes = latitudes(idx);
longitudes = longitudes(idx);

difference_data = [difference_data; ...
                  repmat(150,1,size(difference_data,2));...
                  repmat(-60,1,size(difference_data,2))];
latitudes = [latitudes; -89; -89 ];
longitudes = [longitudes; 0; 0];

% Loop through each comparison to create a map with a boxplot
fonsiz = 24;
scrsz = get(0, 'ScreenSize');

for method_idx = 1:length(method_labels)
    % Extract difference data for the current method
    differences = difference_data(:, method_idx);
    
    % Remove entries with NaN differences
    valid_idx = ~isnan(differences);
    valid_latitudes = latitudes(valid_idx);
    valid_longitudes = longitudes(valid_idx);
    valid_differences = differences(valid_idx);
    
    % Create a new figure with two subplots
    figure('Position', [1 1 scrsz(3) scrsz(4)]);
    
    % Plot the world map in the first subplot
    subplot(1, 6, [1 2 3 4 5]);
    worldmap('World'); % Create a world map
    load coastlines; % Load coastlines for reference
    plotm(coastlat, coastlon, 'k'); % Plot coastlines
    hold on;
    scatterm(valid_latitudes, valid_longitudes, 50, valid_differences, 'filled'); % Color-coded differences
    colormap(custom_colormap); % Change colormap to jet for heatmap
    c = colorbar;
    c.Label.String = 'Difference in Threshold (mm)';
    c.Label.FontSize = fonsiz;
    title(['Difference Map: ', method_labels{method_idx}], 'FontSize', fonsiz, 'Interpreter', 'latex');
    setm(gca, 'FontSize', 12);
    gridm('on'); % Add grid lines
    mlabel('on'); % Add meridian labels
    plabel('on'); % Add parallel labels
    
    % Plot the boxplot in the second subplot
    subplot(1, 6, 6);
    h=boxplot(valid_differences, 'Labels', {method_labels{method_idx}}, 'Colors', 'k','Widths', 0.5);
    
    % Customize outliers to be black
    outliers = findobj(gca, 'Tag', 'Outliers');
    set(outliers, 'MarkerEdgeColor', 'k');  % Set outlier color to black

% Customize boxplot appearance
    set(h, 'LineWidth', 1.5);
    colors = [0.8, 0.8, 0.8]; % Light gray
    boxes = findobj(gca, 'Tag', 'Box');
    for j = 1:length(boxes)
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.5);
    end
    
    % Adjust the size of the boxplot
    % Remove the x-axis and its label
    ax = gca; % Get current axes
    ax.XColor = 'none'; % Hide x-axis line and ticks
    ax.XLabel.String = ''; % Remove x-axis label
    ax.Position = [1.1*ax.Position(1), ax.Position(2)*2.5, ax.Position(3), ax.Position(4)*0.65]; % Adjust [left, bottom, width, height]

    
    grid on;
    set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex' );
    
    % Adjust the figure to remove 20% of the width from the left
    original_position = get(gcf, 'Position'); % Get current figure size and position
    left_trim = original_position(3) * 0.15; % Calculate 20% of the figure width
    % Calculate new height and adjust top and bottom margins
    new_height = original_position(4) * 0.6; % Retain 60% of the original height
    new_bottom = original_position(2) + (original_position(4) * 0.2); % Shift bottom up by 20% of height

    
    
    % Update the figure position to shift it left and reduce the width
    set(gcf, 'Position', [original_position(1) + left_trim, new_bottom, ...
                          original_position(3) - left_trim, new_height]);

    

% % %     
    % Save the combined figure
    saveas(gcf, ['figures\WorldMap_and_Boxplot_' strrep(strrep(method_labels{method_idx}, ' ', '_'), '\%', '') '.png']);
end
close all

%% Define the range of the color bar
color_range = [-50, 150];

% Total number of colors in the colormap
n_colors = 256;

% Position of 0 in the color scale
midpoint = abs(color_range(1)) / (abs(color_range(1)) + color_range(2));

% Number of colors for the negative and positive parts
n_negative = round(midpoint * n_colors);
n_positive = n_colors - n_negative;

% Generate colors for the negative range (from white to blue)
red_negative = linspace(1, 0, n_negative)'; % Red component decreases from 1 to 0
green_negative = linspace(1, 0, n_negative)'; % Green component decreases from 1 to 0
blue_negative = ones(n_negative, 1); % Blue component remains constant at 1
negative_colors = [red_negative, green_negative, blue_negative];

% Generate colors for the positive range (from white to red)
red_positive = ones(n_positive, 1); % Red component remains constant at 1
green_positive = linspace(1, 0, n_positive)'; % Green component decreases from 1 to 0
blue_positive = linspace(1, 0, n_positive)'; % Blue component decreases from 1 to 0
positive_colors = [red_positive, green_positive, blue_positive];

% Combine the negative and positive colors
custom_colormap = [flipud(negative_colors); positive_colors];

%% Differences and pcolor
%Load the station data
load('data\selected_stations_data.mat'); % Contains latitudes, longitudes

% Define threshold data columns for comparison
threshold_data = [data.Langousis_Min1(data.Siglevel==0.05), ...
                  data.Langousis_Min2(data.Siglevel==0.05), ...
                  data.Studentized_Residuals(data.Siglevel==0.01), ...
                  data.Studentized_Residuals(data.Siglevel==0.05), ...
                  data.Anderson_Darling(data.Siglevel==0.01), ...
                  data.Anderson_Darling(data.Siglevel==0.05), ...
                  data.CVM(data.Siglevel==0.01), ...
                  data.CVM(data.Siglevel==0.05)];

% Define method labels
method_labels = {'L2 - L1', ...
                 'SR 1\% - L1', ...
                 'SR 5\% - L1', ...
                 'AD 1\% - L1', ...
                 'AD 5\% - L1', ...
                 'CVM 1\% - L1', ...
                 'CVM 5\% - L1'};

% Extract Langousis Min1 as the baseline
langousis_min1 = threshold_data(:, 1);

% Calculate differences for each method
difference_data = [threshold_data(:, 2) - langousis_min1, ...
                   threshold_data(:, 3) - langousis_min1, ...
                   threshold_data(:, 4) - langousis_min1, ...
                   threshold_data(:, 5) - langousis_min1, ...
                   threshold_data(:, 6) - langousis_min1, ...
                   threshold_data(:, 7) - langousis_min1, ...
                   threshold_data(:, 8) - langousis_min1];

% Filter latitudes and longitudes
latitudes = latitudes(idx);
longitudes = longitudes(idx);

difference_data = [difference_data; ...
                  repmat(150,1,size(difference_data,2));...
                  repmat(-60,1,size(difference_data,2))];
latitudes = [latitudes; -88; -86 ];
longitudes = [longitudes; 0; 0];

% Define grid resolution and maximum interpolation distance
grid_spacing = 2; % Degrees
max_distance = 5; % Maximum allowed distance (degrees)

% Create a global grid
[grid_lon, grid_lat] = meshgrid(-180:grid_spacing:180, -90:grid_spacing:90);

% Loop through each comparison to create a map with a boxplot
fonsiz = 24;
scrsz = get(0, 'ScreenSize');

for method_idx = 1:length(method_labels)
    % Extract difference data for the current method
    differences = difference_data(:, method_idx);
    
    % Remove entries with NaN differences
    valid_idx = ~isnan(differences);
    valid_latitudes = latitudes(valid_idx);
    valid_longitudes = longitudes(valid_idx);
    valid_differences = differences(valid_idx);
    
    % Interpolate data onto the grid
    F = scatteredInterpolant(valid_longitudes, valid_latitudes, valid_differences, 'linear', 'none');
    grid_values = F(grid_lon, grid_lat);
    
    % Calculate distances to nearest stations
    station_coords = [valid_longitudes, valid_latitudes];
    grid_coords = [grid_lon(:), grid_lat(:)];
    distances = pdist2(station_coords, grid_coords);
    min_distances = min(distances, [], 1);
    min_distances = reshape(min_distances, size(grid_lon));
    
    % Mask out regions where distance exceeds max_distance
    grid_values(min_distances > max_distance) = NaN;
    
    % Create a new figure with two subplots
    figure('Position', [1 1 scrsz(3) scrsz(4)]);
    
    % Plot the world map with pcolorm
    subplot(1, 6, [1 2 3 4 5]);
    worldmap('World'); % Create a world map
    load coastlines; % Load coastlines for reference
    plotm(coastlat, coastlon, 'k'); % Plot coastlines
    hold on;
    pcolorm(grid_lat, grid_lon, grid_values); % Interpolated heatmap
    colormap(custom_colormap); % Change colormap to jet for heatmap
    c = colorbar;
    c.Label.String = 'Difference in Threshold (mm)';
    c.Label.FontSize = fonsiz;
    title(['Difference Map: ', method_labels{method_idx}], 'FontSize', fonsiz, 'Interpreter', 'latex');
    setm(gca, 'FontSize', 12);
    gridm('on'); % Add grid lines
    mlabel('on'); % Add meridian labels
    plabel('on'); % Add parallel labels
    
    % Plot the boxplot in the second subplot
    subplot(1, 6, 6);
    h=boxplot(valid_differences, 'Labels', {method_labels{method_idx}}, 'Colors', 'k','Widths', 0.5);
    
    % Customize outliers to be black
    outliers = findobj(gca, 'Tag', 'Outliers');
    set(outliers, 'MarkerEdgeColor', 'k');  % Set outlier color to black

    % Customize boxplot appearance
    set(h, 'LineWidth', 1.5);
    colors = [0.8, 0.8, 0.8]; % Light gray
    boxes = findobj(gca, 'Tag', 'Box');
    for j = 1:length(boxes)
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors, 'FaceAlpha', 0.5);
    end
    
    % Adjust the size of the boxplot
    % Remove the x-axis and its label
    ax = gca; % Get current axes
    ax.XColor = 'none'; % Hide x-axis line and ticks
    ax.XLabel.String = ''; % Remove x-axis label
    ax.Position = [1.1*ax.Position(1), ax.Position(2)*2.5, ax.Position(3), ax.Position(4)*0.65]; % Adjust [left, bottom, width, height]

    
    grid on;
    set(gca, 'FontSize', fonsiz, 'FontName', 'Times', 'TickLabelInterpreter', 'latex' );
    
    % Adjust the figure to remove 20% of the width from the left
    original_position = get(gcf, 'Position'); % Get current figure size and position
    left_trim = original_position(3) * 0.15; % Calculate 20% of the figure width
    % Calculate new height and adjust top and bottom margins
    new_height = original_position(4) * 0.6; % Retain 60% of the original height
    new_bottom = original_position(2) + (original_position(4) * 0.2); % Shift bottom up by 20% of height

    
    
    % Update the figure position to shift it left and reduce the width
    set(gcf, 'Position', [original_position(1) + left_trim, new_bottom, ...
                          original_position(3) - left_trim, new_height]);
    
    % Save the combined figure
    saveas(gcf, ['figures\WorldCotourMap_and_Boxplot_' strrep(strrep(method_labels{method_idx}, ' ', '_'), '\%', '') '.png']);
end

