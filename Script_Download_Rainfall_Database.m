% README FILE FOR DAILY GLOBAL HISTORICAL CLIMATOLOGY NETWORK (GHCN-DAILY) 
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

%% Step 0: URL for the GHCN version file
version_url = 'https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-version.txt';
version_file = fullfile(data_folder, 'ghcnd-version.txt');

% Download the version file
if exist(version_file, 'file') == 0
    websave(version_file, version_url);
end

% Read the version file
fid = fopen(version_file, 'r');
version_info = fgetl(fid);
fclose(fid);

% Display the dataset version
disp(['GHCN-Daily Dataset Version: ', version_info]);

%% Step 1: Download and process the inventory to filter stations with more than 110 years of data
inventory_url = 'https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt';
inventory_file = fullfile(data_folder, 'ghcnd-inventory.txt');

% Check if the file already exists, otherwise download
if exist(inventory_file,'file') == 0
    websave(inventory_file, inventory_url);
    disp('Downloaded inventory file');
end

% Load the inventory file and filter for stations with more than 110 years of data
inventory = readtable(inventory_file, 'Format', '%11s %9f %10f %4s %4d %4d', 'ReadVariableNames', false);
inventory.Properties.VariableNames = {'ID', 'LATITUDE', 'LONGITUDE', 'ELEMENT', 'FIRSTYEAR', 'LASTYEAR'};

% Filter for precipitation data (PRCP) and stations with more than 110 years of records
precipitation_stations = inventory(strcmp(inventory.ELEMENT, 'PRCP') & ...
                                  (inventory.LASTYEAR - inventory.FIRSTYEAR >= 110), :);

%% Step 2: Download .dly Files and Process Precipitation Data
base_url = 'https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/';

% Loop over all station IDs and download their corresponding data
for i = 1:height(precipitation_stations)
    station_id = precipitation_stations.ID{i};
    station_file = [station_id, '.dly'];
    station_url = [base_url, station_file];
    
    % Define paths for the station file and output CSV
    station_path = fullfile(data_folder, station_file);
    output_file = fullfile(data_folder, [station_id, '_precipitation.csv']);
    
    % Check if the station file already exists locally
    if exist(output_file, 'file')
        disp(['File already exists locally: ', output_file]);
    else
        % Download the .dly file
        try
            websave(station_path, station_url);
            disp(['Downloaded file: ', station_file]);
        catch
            disp(['Failed to download data for station ', station_id]);
            continue;
        end

        % Process the .dly file to extract date and precipitation records
        fid = fopen(station_path);
        precip_data = {};
        while ~feof(fid)
            line = fgetl(fid);
            if strcmp(line(18:21), 'PRCP')  % Filter for precipitation data
                year = line(12:15);
                month = line(16:17);
                for day = 1:31
                    value = str2double(line(22+(day-1)*8:26+(day-1)*8));
                    if value ~= -9999
                        date = datenum([year, '-', month, '-', num2str(day, '%02d')], 'yyyy-mm-dd');
                        precip_data = [precip_data; {datestr(date, 'yyyy-mm-dd'), value / 10}];  % Precipitation in mm
                    end
                end
            end
        end
        fclose(fid);
        
        % Save the extracted data to a file
        precip_table = cell2table(precip_data, 'VariableNames', {'Date', 'Precipitation_mm'});
        writetable(precip_table, output_file);
        disp(['Processed data for station ', station_id, ' and saved to ', output_file]);
        
        % Delete the .dly file after processing
        if exist(station_path, 'file') == 2
            delete(station_path);
            disp(['Deleted file: ', station_path]);
        else
            disp(['File not found: ', station_path]);
        end
    end
end

%% Step 3: Calculate data availability and Save to filelistnames.txt
output_filelist = fullfile(data_folder, 'filelistnames.txt');
if exist(output_filelist, 'file')
    disp(['File already exists locally: ', output_filelist]);
else
    fileID = fopen(output_filelist, 'w');
    fprintf(fileID, 'Filename, Percentage Available, Total Days, Available Rows\n');
    
    % Process each CSV file
    csv_files = dir(fullfile(data_folder, '*_precipitation.csv'));
    for i = 1:length(csv_files)
        filename = fullfile(data_folder, csv_files(i).name);
        
        % Load the data
        data = readtable(filename);
        
        % Calculate total days and available rows
        first_date = min(datenum(data.Date, 'yyyy-mm-dd'));
        last_date = max(datenum(data.Date, 'yyyy-mm-dd'));
        total_days = last_date - first_date + 1;
        available_rows = height(data);
        percentage_available = (available_rows / total_days) * 100;
        
        % Append the result to the file
        fprintf(fileID, '%s, %.2f, %d, %d\n', csv_files(i).name, percentage_available, total_days, available_rows);
    end
    fclose(fileID);
    disp('File list with data availability saved to filelistnames.txt');
end

%% Step 4: Plot a histogram of 'Percentage Available'
data = readtable(output_filelist);
% Create the histogram with gray fill, black edges, and LaTeX fonts
fonsiz = 18;
scrsz = get(0, 'ScreenSize');
figure('Position', [1 1 scrsz(3) scrsz(4)]);
ax_ = newplot;
histogram(data.PercentageAvailable, 'BinWidth', 5, 'EdgeColor', 'black', 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.7);

title('Histogram of Files per Percentage of data available ', 'Interpreter', 'latex');
xlabel('Percentage of data available (\%)', 'Interpreter', 'latex');
ylabel('Frequency of files', 'Interpreter', 'latex');
grid on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fonsiz, 'TickLabelInterpreter', 'latex');

% Save the histogram in the figures folder
saveas(gcf, fullfile(figures_folder, 'HistogramQualityRecords.png'));

%% Step 5: Generate Threshold Count Plot
thresholds = 0.99:-0.01:0.05;
counts = arrayfun(@(x) sum(data.PercentageAvailable >= x * 100), thresholds);

% Plot threshold count plot
figure('Position', [1 1 scrsz(3) scrsz(4)]);
plot(thresholds * 100, counts, '-', 'Color', 'k', 'LineWidth', 1.5);  % Black line for the plot

title('Count of Files per Percentage of data available', 'Interpreter', 'latex');
xlabel('Percentage of data available (\%)', 'Interpreter', 'latex');
ylabel('Count', 'Interpreter', 'latex');
set(gca, 'XTick', 0:5:100, 'YTick', 0:250:max(counts) + 250);
grid on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fonsiz, 'TickLabelInterpreter', 'latex');

% Mark the 90% threshold point
threshold_90 = 90;
index_90 = find(thresholds * 100 == threshold_90);
count_90 = counts(index_90);
hold on;
plot(threshold_90, count_90, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(threshold_90 - 30, count_90, sprintf('90%% (%d files)', count_90), ...
    'FontSize', fonsiz, 'FontName', 'Helvetica', 'Interpreter', 'latex', 'Color', 'r');

% Save the threshold count plot in the figures folder
saveas(gcf, fullfile(figures_folder, 'ECDFQualityRecords.png'));

