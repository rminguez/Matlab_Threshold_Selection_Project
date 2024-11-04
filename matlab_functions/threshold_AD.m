function threshold = threshold_AD(data, alpha)
    % threshold_AD computes the optimal threshold using the Anderson-Darling test.
    % The function fits a Generalized Pareto Distribution (GPD) to the data's excesses 
    % above different threshold candidates and uses the Anderson-Darling (AD) test to select the best threshold.
    %
    % Inputs:
    % - data: A vector of data points to analyze.
    % - alpha: (Optional) The significance level for the Anderson-Darling test (default is 0.05).
    %
    % Output:
    % - threshold: The optimal threshold found based on the AD test.
    
    % Set the default alpha value if not provided
    if nargin < 2
        alpha = 0.05;
    end
    
    % Sort the data and remove duplicate values to define unique threshold candidates
    candidates = sort(unique(data(:)));  % Ensure data is sorted and free of duplicates
    
    % Initialize the p-value and index for the threshold search
    p_value = 0;  % Initialize p-value for the Anderson-Darling test
    i = 1;  % Start with the first candidate threshold
    
    % Loop through the threshold candidates until p-value exceeds alpha
    while p_value < alpha
        % Set the current threshold to the i-th candidate
        u = candidates(i);
        
        % Calculate excesses: values in the data greater than the threshold u
        excesses = data(data > u) - u;
        
        % If there are fewer than 10 excesses, stop the search (insufficient data for GPD fitting)
        if length(excesses) < 10
            break;
        end
        
        % Fit a Generalized Pareto Distribution (GPD) to the excesses using maximum likelihood estimation (MLE)
        paramEsts = mle(excesses, 'distribution', 'gp');
        scale = paramEsts(1);  % Scale parameter (sigma)
        shape = paramEsts(2);  % Shape parameter (xi)
        
        % Anderson-Darling test: Calculate the CDF of the fitted GPD
        PGPD = @(q) gpcdf(q, shape, scale, u);  % Define CDF of the fitted GPD
        
        % Compute the probabilities for each excess
        prob = PGPD(excesses);
        % Remove probabilities too close to 0 or 1 to avoid numerical issues in the AD test
        prob = prob(prob > 1e-8 & prob < 1 - 1e-8);
        n = length(prob);  % Number of valid probabilities
        k = (1:n)';  % Indices for each data point
        
        % Transform the probabilities into normal quantiles
        qnor = norminv(sort(prob), 0, 1);  % Normal inverse CDF (Z-scores)
        % Normalize the quantiles to calculate the adjusted probabilities (pnor)
        pnor = normcdf((qnor - mean(qnor)) / std(qnor));  % Standard normal CDF
        
        % Compute Anderson-Darling statistic (A)
        A = (-n - sum((2 * k - 1) .* log(pnor) + (2 * n + 1 - 2 * k) .* log(1 - pnor)) / n);
        A = round((1 + 0.75 / n + 2.25 / n^2) * A, 4);  % Adjust the AD statistic for small sample sizes
        
        % Calculate the p-value based on the Anderson-Darling statistic
        if A < 0.2
            p_value = 1 - exp(-13.436 + 101.14 * A - 223.73 * A^2);
        elseif A >= 0.2 && A < 0.34
            p_value = 1 - exp(-8.318 + 42.796 * A - 59.938 * A^2);
        elseif A >= 0.34 && A < 0.6
            p_value = exp(0.9177 - 4.279 * A - 1.38 * A^2);
        elseif A >= 0.6 && A < 10
            p_value = exp(1.2937 - 5.709 * A + 0.0186 * A^2);
        else
            p_value = 7.37e-10;  % If A >= 10, set a very small p-value
        end
        
        % Move to the next threshold candidate
        i = i + 1;
        
        % Stop the loop if all candidates have been tested
        if i > length(candidates)
            break;
        end
    end
    
    % Return the last threshold found before the p-value exceeded alpha
    threshold = u;
end
