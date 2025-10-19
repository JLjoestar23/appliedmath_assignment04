% GLOBAL_TRUNC_ERROR_ANALYSIS_MODIFIED
%   Plots the global truncation error versus both step size (h)
%   and total number of function evaluations on log–log scales
%   for one or more numerical integration methods. Also estimates
%   the order of accuracy from the slope of a linear fit in log–log space.
%
% INPUTS:
%   method_list     - Cell array of method names 
%                     (e.g., {'ForwardEuler','MidpointMethod'})
%                     Each must correspond to a valid input for
%                     calc_global_trunc_error().
%
%   rate_func       - Function handle describing system dynamics:
%                       dXdt = rate_func(t, X)
%
%   analytical_soln - Function handle returning the true solution X(t),
%                     used to compute the global truncation error.
%
%   tspan           - Two-element vector [t_start, t_end] defining the
%                     integration interval.
%
%   num_trials      - Number of step sizes (h) to test. Step sizes are
%                     logarithmically spaced between 1e-6 and 10^1.
%
% OUTPUTS:
%   Two log–log plots:
%     1. Global Truncation Error vs. Step Size (h)
%     2. Global Truncation Error vs. Number of Function Calls

function global_trunc_error_analysis_modified(method_list, rate_func, analytical_soln, tspan, num_trials)

    % generate log-spaced step sizes
    h_list = logspace(-4, 1, num_trials);
    h_avg_list = zeros(length(method_list), num_trials);

    % initialize storage arrays
    global_errors = zeros(length(method_list), num_trials);
    num_evals_all = zeros(length(method_list), num_trials);

    % plot style presets
    error_data_presets = {'b.', 'r.', 'm.', 'k.'};

    % main loop through each integration method
    for i = 1:length(method_list)
        method_name = method_list{i};

        % compute global truncation error for each step size
        for j = 1:num_trials
            [global_errors(i, j), h_avg_list(i, j), num_evals_all(i, j)] = ...
                calc_global_trunc_error(method_name, rate_func, analytical_soln, tspan, h_list(j));
        end

        % plot: global truncation error vs. step size
        figure(1);
        loglog(h_avg_list(i,:), global_errors(i,:), error_data_presets{i}, ...
               'MarkerSize', 5, 'DisplayName', method_name);
        hold on;

        % second plot: global truncation error vs. number of evals
        figure(2);
        loglog(num_evals_all(i,:), global_errors(i,:), error_data_presets{i}, ...
               'MarkerSize', 5, 'DisplayName', sprintf('%s Data', method_name));
        hold on;
    end

    % plot 1: error vs. step size
    figure(1);
    title('Global Truncation Error vs. Step Size');
    xlabel('Timestep (s)');
    ylabel('Global Truncation Error');
    legend('Location', 'northwest');
    grid on;
    axis tight;
    hold off;

    % plot 2: error vs. number of function evaluations
    figure(2);
    title('Global Truncation Error vs. Number of Function Calls');
    xlabel('Number of Function Calls');
    ylabel('Global Truncation Error');
    legend('Location', 'northeast');
    grid on;
    axis tight;
    hold off;
end
