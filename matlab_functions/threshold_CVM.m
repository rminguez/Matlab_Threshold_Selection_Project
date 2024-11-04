function threshold = threshold_CVM(data, alpha)
    % threshold_CVM computes the optimal threshold based on the Cramer-Von Mises (CVM) statistic.
    % The function fits a Generalized Pareto Distribution (GPD) to the data's excesses 
    % above different threshold candidates and uses the CVM test to select the best threshold.
    %
    % Inputs:
    % - data: A vector of data points to analyze.
    % - alpha: (Optional) The significance level for the CVM test (default is 0.05).
    %
    % Output:
    % - threshold: The optimal threshold found based on the CVM test.

    % Set the default alpha value if not provided
    if nargin < 2
        alpha = 0.05;
    end
    
    % Sort the data and retain only unique values as potential threshold candidates
    candidates = sort(unique(data(:)));  % Remove duplicate values and sort in ascending order
    
    % Initialize the p-value and index
    p_value = 0;  % Initialize p-value for the CVM test
    i = 1;  % Start with the first candidate threshold
    
    % Continue testing thresholds until p-value exceeds the significance level (alpha)
    while p_value < alpha
        % Define the current threshold
        u = candidates(i);
        
        % Compute the excesses (data points above the current threshold)
        excesses = data(data > u) - u;
        
        % If there are fewer than 10 excesses, stop the search (insufficient data)
        if length(excesses) < 10
            break;
        end
        
        % Fit a Generalized Pareto Distribution (GPD) to the excesses using maximum likelihood estimation (MLE)
        paramEsts = mle(excesses, 'distribution', 'gp');
        scale = paramEsts(1);  % Scale parameter (sigma)
        shape = paramEsts(2);  % Shape parameter (xi)
        
        % Define the CDF of the fitted GPD
        PGPD = @(q) gpcdf(q, shape, scale, u);
        
        % Calculate the cumulative probabilities (CDF) of the excesses
        prob = PGPD(excesses);
        % Remove probabilities very close to 0 and 1 to avoid numerical issues
        prob = prob(prob > 1e-8 & prob < 1 - 1e-8);
        n = length(prob);  % Number of valid probabilities
        k = (1:n)';  % Indices for each data point
        
        % Transform the probabilities to the normal space using the inverse normal CDF
        qnor = norminv(sort(prob));
        % Recalculate the probabilities based on the standard normal distribution
        pnor = normcdf((qnor - mean(qnor)) / std(qnor));
        
        % Compute the Cramer-Von Mises (CVM) statistic
        w = round((sum((pnor - (2 * k - 1) / (2 * n)).^2) + 1 / (12 * n)) * (1 + 0.5 / n), 4);
        
        % Calculate the p-value based on the CVM statistic
        if w < 0.0275
            p_value = 1 - exp(-13.953 + 775.5 * w - 12542.61 * w^2);
        elseif w >= 0.0275 && w < 0.051
            p_value = 1 - exp(-5.903 + 179.546 * w - 1515.29 * w^2);
        elseif w >= 0.051 && w < 0.092
            p_value = exp(0.886 - 31.62 * w + 10.897 * w^2);
        elseif w >= 0.092 && w < 1.1
            p_value = exp(1.111 - 34.242 * w + 12.832 * w^2);
        elseif w >= 1.1
            p_value = 7.37e-10;
        end
        
        % Move to the next threshold candidate
        i = i + 1;
        
        % If all candidates have been tested, exit the loop
        if i > length(candidates)
            break;
        end
    end
    
    % Return the optimal threshold found
    threshold = u;
end
