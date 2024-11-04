function local_minima = threshold_MSE(pks_unicos_valid, excedencias_mean_valid, excedencias_weight_valid, n0, p, plot_flag, filename, display_flag)
    % threshold_Langousis computes Langousis statistics and finds candidate thresholds
    % based on local minima of a spline fit. Optionally, plots the graph if plot_flag is set.
    %
    % Inputs:
    % - pks_unicos_valid: vector of unique peaks (potential thresholds)
    % - excedencias_mean_valid: vector of exceedance means
    % - excedencias_weight_valid: vector of exceedance weights
    % - n0: (optional) number of points to skip from the end (default = 10)
    % - p: (optional) smoothing parameter for cubic spline fitting (default = 0.9)
    % - plot_flag: (optional) boolean flag, true to plot the graph, false to skip plotting (default = false)
    % - filename: (optional) path and name for making graphs
    % - display_flag: (optional) boolean flag, true to display messages, false otherwise
    %
    % Output:
    % - local_minima: candidate thresholds (local minima of the spline)

    % Set default values if not provided
    if nargin < 4
        n0 = 10;
    end
    if nargin < 5 || isempty(p)
        p = 0.9;
    end
    if nargin < 6 || isempty(plot_flag)
        plot_flag = false;
    end
    if nargin < 7 || isempty(filename)
        filename = [];
    end
    if nargin < 8 || isempty(display_flag)
        display_flag = false;
    end

    % Initialize Langousis and chi-square test values
    Langousis = zeros(length(pks_unicos_valid) - n0, 1);
    
    % Compute Langousis and chi-square test values
    for i = 1:length(pks_unicos_valid) - n0
        [~, fobj] = RWLSfit(pks_unicos_valid(i:end), excedencias_mean_valid(i:end), excedencias_weight_valid(i:end));
        Langousis(i) = fobj / length(pks_unicos_valid(i:end));
    end

    % Define x and y for spline fitting
    x = pks_unicos_valid(1:length(pks_unicos_valid) - n0);   % Thresholds
    y = log(Langousis);  % Logarithm of Langousis values

    % Step 1: Fit a cubic spline to the data
    spline_fit = csaps(x, y, p);  % Smoothing parameter p

    % Step 2: Differentiate the spline to get its first derivative
    spline_derivative = fnder(spline_fit);

    % Step 3: Find the roots of the derivative (where it is zero)
    extrema = fnzeros(spline_derivative, [min(x), max(x)]);

    % Step 4: Compute second derivative to classify minima
    spline_second_derivative = fnder(spline_derivative);

    % Evaluate the second derivative at the points where the first derivative is zero
    extrema_vals = extrema(1, :);  % Extract x values of extrema
    second_derivative_vals = fnval(spline_second_derivative, extrema_vals);

    % Local minima occur where the second derivative is positive
    local_minima = extrema_vals(second_derivative_vals > 0);

    % Display local minima
    if display_flag
        disp('Local Minima at:');
        disp(local_minima);
    end

    % Step 5 (Optional): Plot if plot_flag is true
    if plot_flag
        
        fonsiz = 18;
        scrsz = get(0, 'ScreenSize');

        % Create a new figure with screen size
        figure('Position', [1 1 scrsz(3) scrsz(4)]);
        ax_ = newplot;
        legh_ = [];
        legt_ = {};

        % Plot Langousis values (MSE data) in semilogarithmic scale
        h_ = semilogy(x, Langousis, 'o', 'LineWidth', 2, 'MarkerSize', 5,'color',[0.5 0.5 0.5]);  % Plot Langousis values
        legh_(end + 1) = h_;
        legt_{end + 1} = 'MSE data';
        hold on;

        % Plot the spline fit
        spline_plot_data = fnplt(spline_fit);  % Get spline points (X and Y)
        h_ = plot(spline_plot_data(1, :), exp(spline_plot_data(2, :)), 'k-', 'LineWidth', 2);  % Plot the spline (use exp to reverse log scale)
        legh_(end + 1) = h_;
        legt_{end + 1} = 'Spline Fit';

        % Plot the local minima
        h_ = plot(local_minima, exp(fnval(spline_fit, local_minima)), 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 10);  % Plot local minima
        legh_(end + 1) = h_;
        legt_{end + 1} = 'Local Minima';

        % Add gray dashed vertical lines from the local minima to the bottom of the graph
        for i = 1:length(local_minima)
            plot([local_minima(i) local_minima(i)], [min(exp(spline_plot_data(2, :))) max(Langousis)], 'k--', 'Color', [0.5 0.5 0.5]);

            % Add rotated text at each local minima
            hh = text(local_minima(i) + 0.02 * range(x), min(exp(spline_plot_data(2, :))) * 1.1, ...
                ['Threshold = ' num2str(local_minima(i))], 'Rotation', 90, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
        end
        
        % Grid and labels
        grid on;
        xlabel('Thresholds (u)', 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
        ylabel('Mean Leasts Squares Error', 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
        title('Mean Leasts Squares Error Values and Spline Fit', 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');

        % Add the legend
        legend(legh_, legt_, 'Location', 'NorthEast', 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');

        % Format the axes
        set(gca, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'TickLabelInterpreter', 'latex');
        set(gcf, 'PaperPositionMode', 'auto');

        hold off;
        
        if ~isempty(filename)
            saveas(gcf, [filename 'Langousis'], 'png');
            saveas(gcf, [filename 'Langousis'], 'epsc');
        end

    end
end
