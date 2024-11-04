function [threshold,beta,fobj,r] = threshold_studentized_residuals(pks_unicos_valid, excedencias_mean_valid, excedencias_weight_valid, siglevel, plot_flag, filename, display_flag)
    % threshold_studentized_residuals computes the optimal threshold based on Chi-squared
    % and studentized residuals. Optionally plots the results if plot_flag is true and 
    % displays messages if display_flag is true.
    %
    % Inputs:
    % - pks_unicos_valid: vector of unique peaks (potential thresholds)
    % - excedencias_mean_valid: vector of exceedance means
    % - excedencias_weight_valid: vector of exceedance weights
    % - siglevel: (optional) significance level for Chi-squared test (default 0.05)
    % - plot_flag: (optional) boolean flag, true to plot the graphs, false otherwise
    % - filename: (optional) path and name for making graphs
    % - display_flag: (optional) boolean flag, true to display messages, false otherwise
    %
    % Output:
    % - threshold: the optimal threshold found
    % - beta: optimal regression coefficients
    % - fobj: optimal objective function (weighted leats squares)
    % - r: optimal residuals

    % Default values for optional arguments
    if nargin < 4 || isempty(siglevel)
        siglevel = 0.05;
    end
    if nargin < 5 || isempty(plot_flag)
        plot_flag = false;
    end
    if nargin < 6 || isempty(filename)
        filename = [];
    end
    if nargin < 7 || isempty(display_flag)
        display_flag = false;
    end

    stop_search = 0;
    it = 1;
    threshold = pks_unicos_valid(1);  % Initial threshold

    while ~stop_search && it <= 10

        % Find the current threshold in the pks_unicos_valid array
        pos = find(pks_unicos_valid == threshold);
        u_values = pks_unicos_valid(pos:end);  % Thresholds starting from the current one
        e_values = excedencias_mean_valid(pos:end);  % Exceedances
        w_values = excedencias_weight_valid(pos:end);  % Weights

        % Perform the RWLS fitting and calculate studentized residuals
        [beta, fobj, r, rN] = RWLSfit(u_values, e_values, w_values);

        % Plot residuals if requested
        if plot_flag
            fonsiz = 18;
            scrsz = get(0, 'ScreenSize');
            figure('Position', [1 1 scrsz(3) scrsz(4)]);
            ax_ = newplot;
            legh_ = [];
            legt_ = {};
            h_ = plot(u_values, rN, 'k', 'LineWidth', 2);
            legh_(end + 1) = h_;
            legt_{end + 1} = ['Internally studentized residuals'];
            hold on;
            grid on;
            hh = xlabel('Threshold $u$ (mm/d)');
            set(hh, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
            hh = ylabel('$r$');
            set(hh, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
            hh = text(5 + threshold, min(rN) + 0.1 * (max(rN) - min(rN)), ['Min threshold = ' num2str(threshold)]);
            set(hh, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
            hold off;
            leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
            h_ = legend(ax_, legh_, legt_, leginfo_{:});  % create legend
            set(h_, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
            set(gca, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'TickLabelInterpreter', 'latex');
            set(gcf, 'PaperPositionMode', 'auto');
            if ~isempty(filename)
                saveas(gcf, [filename 'StudenRes' num2str(it)], 'png');
                saveas(gcf, [filename 'StudenRes' num2str(it)], 'epsc');
            end
        end

        % Check stopping criteria: Chi-squared test and studentized residuals
        if fobj > chi2inv(1 - siglevel, length(u_values) - 2) || abs(rN(1)) > norminv(1 - siglevel / 2, 0, 1)
            if display_flag
                if fobj > chi2inv(1 - siglevel, length(u_values) - 2)
                    disp('Chi-squared test detects anomalies');
                end
                if abs(rN(1)) > norminv(1 - siglevel / 2, 0, 1)
                    disp('The maximum studentized residual of the first record detects anomalies');
                end
            end
            thresholdsearch = 1;
        else
            thresholdsearch = 0;
            stop_search = 1;  % If criteria met, stop the loop
        end

        % If anomalies detected, perform threshold search
        if thresholdsearch
            if display_flag
                disp(['Maximum sensitivity  = ' num2str(max(abs(rN))) ' and thus the optimal threshold seems to be on the right side of the minimum sample value, looking for the location']);
            end
            [~, threshold] = threshold_search(u_values, rN, w_values, plot_flag, [filename 'thresholdlocation' num2str(it)]);
            if display_flag
                disp(['New threshold found: ' num2str(threshold)]);
            end
        end

        it = it + 1;
    end

    % Subfunction: threshold_search
    function [fitresult, threshold] = threshold_search(u_data, e_data, W_data, ploteat, filename)
        % Renamed function variables to avoid conflict with outer function
        %
        % Inputs:
        % - u_data: threshold values
        % - e_data: exceedances
        % - W_data: weights
        % - ploteat: flag for plotting
        % - filename: file name to save plots
        % Outputs:
        % - fitresult: fit object representing the fit
        % - threshold: the threshold value determined from the fit

        if nargin < 3 || isempty(W_data)
            W_data = ones(size(u_data));
        end
        if nargin < 4 || isempty(ploteat)
            ploteat = 1;
        end
        if nargin < 5 || isempty(filename)
            filename = [];
        end

        e_data = e_data(:);
        u_data = u_data(:);
        [u_data, ord] = sort(u_data);
        e_data = e_data(ord);
        W_data = W_data(ord);

        % Fit: Smoothing spline
        SmoothingParam = fminbnd(@(x) (smoothingspline(u_data, e_data, W_data, x) - 0.9)^2, 0.5, 0.99);
        [~, fitresult, ~] = smoothingspline(u_data, e_data, W_data, SmoothingParam);

        % Find the first zero from the left
        uc = linspace(u_data(1), u_data(end), 1000)';
        ec = fitresult(uc);
        currentsign = sign(ec(1));
        zeroloc = [0 0];
        cont = 1;
        for i = 2:length(ec)
            if currentsign ~= sign(ec(i))
                zeroloc(cont) = (uc(i) + uc(i-1)) / 2;
                cont = cont + 1;
                currentsign = -currentsign;
            end
            if cont > 2
                break;
            end
        end

        pos1 = find(u_data >= zeroloc(1) & u_data <= zeroloc(2));
        [mini, posi] = max(abs(e_data(pos1)));
        posi = pos1(1) + posi - 1;
        threshold = u_data(posi);
        mini = e_data(posi);

        if ploteat
            % Plot fit with data.
            fonsiz = 18;
            scrsz = get(0, 'ScreenSize');
            figure('Position', [1 1 scrsz(3) scrsz(4)]);
            ax_ = newplot;
            legh_ = [];
            legt_ = {};
            h_ = plot(u_data, e_data, '.k', 'MarkerSize', 10);
            legh_(end + 1) = h_;
            legt_{end + 1} = ['Data'];
            grid on;
            hold on;
            plot(threshold * ones(1, 100), linspace(min(e_data), max(e_data), 100), '--', 'color', [0.5 0.5 0.5], 'LineWidth', 2);
            h_ = plot(threshold, mini, 'ok', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'w');
            legh_(end + 1) = h_;
            legt_{end + 1} = ['Local optimum'];
            hh = xlabel('Threshold $u$');
            set(hh, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
            hh = ylabel('$r^N$');
            set(hh, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'Interpreter', 'latex');
            legend(ax_, legh_, legt_, 'Orientation', 'vertical', 'Location', 'NorthEast');
            set(gca, 'FontName', 'Montserrat', 'FontSize', fonsiz, 'TickLabelInterpreter', 'latex');
            set(gcf, 'PaperPositionMode', 'auto');
            if ~isempty(filename)
                saveas(gcf, [filename '.png']);
                saveas(gcf, [filename], 'epsc');
            end
        end

        function [qualityparam, fitresult, gof] = smoothingspline(x_data, y_data, w_data, SmoothingParam)
            % Set up fittype and options for smoothing spline
            ft = fittype('smoothingspline');
            opts = fitoptions('Method', 'SmoothingSpline');
            opts.Normalize = 'on';
            opts.SmoothingParam = SmoothingParam;

            % Fit model to data
            [fitresult, gof] = fit(x_data, y_data, ft, opts);
            qualityparam = gof.rsquare;
        end
    end
end
