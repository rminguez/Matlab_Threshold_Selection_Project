% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %   SCRIPT para optimizar la selección del umbral Registro RCP_ASN00021043
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
actualpath=pwd;

%   Añado las carpetas con los paths
addpath([actualpath '\matlab_functions'])

%%  Graficado de los datos
% %
graficos = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   REGIMEN EXTREMAL ANUAL DE LAS SERIES
%%%%%%%%
if exist(['figures\PaperThreshold'],'dir')~=7,
    mkdir(['figures\PaperThreshold'])
end

ejemplo='PRCP_ASN00021043';

%   Si existe hago la lectura del fichero
if strcmp(ejemplo,'PRCP_ASN00021043'),
    %exist([actualpath '\data\NOAA-NCDCprecipitations\PRCP_ASN00021043.csv'],'file')==2,
    
    data=readtable([actualpath '\data\PRCP_ASN00021043.csv'],'TreatAsEmpty',{'NA'},'format','%s%f');
    % Antes de hacer nada vamos a eliminar los registros nulos o NaN
    % Verificar la segunda columna (precipitación)
    rowsToDelete = isnan(data{:, 2});

    % Eliminar las filas que cumplen la condición
    data(rowsToDelete, :) = [];
    
        %
    pluviometros.nombre = 'PRCP_ASN00021043';
    pluviometros.datenum =datenum(data{:,1});
    pluviometros.fechas =data{:,1};
    pluviometros.data =data{:,2}/10;
    
    % Calcular la distancia mínima entre registros consecutivos en el tiempo
    % Esto se asume que es la distancia en días
    min_time_diff = min(diff(pluviometros.datenum));
    
    % Definir la distancia mínima entre picos en unidades temporales
    % Por ejemplo, la distancia mínima será el doble de la distancia mínima entre registros
    min_peak_distance = 2 * min_time_diff;
    
    %   La ejecucion con un n la guardo en el fichero
    %   workspaceMC_POT_AM_n.mat
% % %     load(['PRCP_ASN00021043.mat'])
    threshold0 = 2.0;
end

%   Defino un threshold inicial
threshold=0.0; 
n0=10;
siglevel=0.01;
%   Extraigo todos los thresholds válidis en orden y sus correspondientes
%   medias de excedencias y pesos
[pks_unicos_valid,excedencias_mean_valid,excedencias_weight_valid, pks, locs] = threshold_peak_extraction(pluviometros.data,threshold,n0,min_peak_distance);

% Opcional: Gráfico de los residuos estudentizados
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

%% Anderson_darling method

threshold_val_AD = threshold_AD(pluviometros.data, siglevel);
threshold_val_AD = threshold_AD(pks, siglevel);

%% Crame-Von Misses
threshold_val_CVM = threshold_CVM(pluviometros.data, siglevel);
threshold_val_CVM = threshold_CVM(pks, siglevel);

[threshold_val_SR,threshold_val_MSE,threshold_val_AD,threshold_val_CVM]


