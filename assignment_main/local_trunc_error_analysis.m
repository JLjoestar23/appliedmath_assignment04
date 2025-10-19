% LOCAL_TRUNC_ERROR_ANALYSIS
%   Plots the local truncation error versus step size (h) on a log–log
%   scale for one or more numerical integration methods, and estimates
%   each method’s order of accuracy from the slope of a line fit.
%
% INPUTS:
%   method_list     - Cell array of method names (e.g., {'ForwardEuler','MidpointMethod'})
%                     Each method name must correspond to a valid input
%                     for the function calc_local_trunc_error().
%
%   rate_func       - Function handle describing the system dynamics:
%                       dXdt = rate_func(t, X)
%
%   analytical_soln - Function handle returning the true solution X(t),
%                     used to compute the truncation error.
%
%   t_ref           - Reference time at which the local truncation error
%                     is evaluated.
%
%   num_trials      - Number of step sizes (h) to test. Step sizes are
%                     logarithmically spaced between 1e-6 and 10^1.
% OUTPUTS:
%   loglog plot representing Global Truncation Error vs. Step Size (h)

function local_trunc_error_analysis(BT_struct_list, rate_func, analytical_soln, t_ref, num_trials)
    
    % initialize step sizes (log-spaced)
    h_list = logspace(3, 7, num_trials); % seconds
    
    % initialize storage
    local_truncation_errors = zeros(length(BT_struct_list), num_trials);
    local_h_diff = zeros(1, length(BT_struct_list));
    
    % plot style presets
    error_data_presets = {'r.', 'b.', 'm.'};
    line_fit_presets = {'g--', 'c--', 'y--'};

    % loop through methods
    for i = 1:length(BT_struct_list)
        BT_struct = BT_struct_list{i};
        
        % compute local truncation errors for this method
        for j = 1:num_trials
            [local_truncation_errors(i, j), local_h_diff(j)] = ...
                calc_local_trunc_error(rate_func, analytical_soln, t_ref, h_list(j), BT_struct);
        end

        % perform linear regression in log–log space (middle 60%)
        start_i = floor(0.3*num_trials);
        end_i   = floor(0.6*num_trials);

        log_h_list = log10(h_list(start_i:end_i));
        log_local_trunc_errors = log10(local_truncation_errors(i, start_i:end_i));

        p = polyfit(log_h_list, log_local_trunc_errors, 1); % [slope, intercept]
        slope = p(1);
        fit_vals = polyval(p, log_h_list);

        % plot actual data
        loglog(h_list, local_truncation_errors(i,:), error_data_presets{i}, ...
               'MarkerSize', 5, 'DisplayName', BT_struct.name);
        hold on;

        % plot linear fit
        loglog(10.^log_h_list, 10.^fit_vals, line_fit_presets{i}, ...
               'LineWidth', 2, ...
               'DisplayName', sprintf('Fit (slope = %.2f)', slope));
    end
    
    % plot local timestep difference
    loglog(h_list, local_h_diff, 'k-', 'LineWidth', 2, 'DisplayName', 'Local Timestep Diff')

    % labels and aesthetics
    title('Local Truncation Error vs. Step Size');
    xlabel('Timestep (s)');
    ylabel('Local Truncation Error');
    legend('Location', 'northwest');
    grid on;
    axis tight;
    hold off;
end