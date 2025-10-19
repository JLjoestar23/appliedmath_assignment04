% LOCAL_TRUNC_ERROR_ANALYSIS_MODIFIED
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
%   loglog plot representing local truncation error vs. step size (h)

function local_trunc_error_analysis_modified(method_list, rate_func, analytical_soln, t_ref, num_trials)
    
    % initialize step sizes (log-spaced)
    h_list = logspace(-6, 1, num_trials);
    
    % initialize storage
    local_truncation_errors = zeros(length(method_list), num_trials);
    local_h_diff = zeros(1, length(num_trials));
    
    % plot style presets
    error_data_presets = {'r.', 'b.', 'm.', 'k.'};

    % loop through methods
    for i = 1:length(method_list)
        method_name = method_list{i};
        
        % compute local truncation errors for this method
        for j = 1:num_trials
            [local_truncation_errors(i, j), local_h_diff(j)] = ...
                calc_local_trunc_error(method_name, rate_func, analytical_soln, t_ref, h_list(j));
        end

        % plot actual data
        loglog(h_list, local_truncation_errors(i,:), error_data_presets{i}, ...
               'MarkerSize', 5, 'DisplayName', method_name);
        hold on;
    end

    % labels and aesthetics
    title('Local Truncation Error vs. Step Size');
    xlabel('Timestep (s)');
    ylabel('Local Truncation Error');
    legend('Location', 'northwest');
    grid on;
    axis tight;
    hold off;
end