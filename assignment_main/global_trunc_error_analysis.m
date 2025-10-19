% GLOBAL_TRUNC_ERROR_ANALYSIS
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

function global_trunc_error_analysis(BT_list, rate_func, analytical_soln, tspan, num_trials)

    % generate log-spaced step sizes
    h_list = logspace(3, 7, num_trials);
    h_avg_list = zeros(length(BT_list), num_trials);

    % initialize storage arrays
    global_errors = zeros(length(BT_list), num_trials);
    num_evals_all = zeros(length(BT_list), num_trials);

    % plot style presets
    error_data_presets = {'b.', 'r.', 'm.'};
    line_fit_presets = {'y--', 'g--', 'c--'};

    % main loop through each integration method
    for i = 1:length(BT_list)
        BT = BT_list{i};

        % compute global truncation error for each step size
        for j = 1:num_trials
            [global_errors(i, j), h_avg_list(i, j), num_evals_all(i, j)] = ...
                calc_global_trunc_error(BT, rate_func, analytical_soln, tspan, h_list(j));
        end

        % choose central range for line fitting (middle 60%)
        start_i = floor(0.2*num_trials);
        end_i   = floor(0.6*num_trials);

        % log-transformed data for regression
        log_h = log10(h_avg_list(i, start_i:end_i));
        log_err = log10(global_errors(i, start_i:end_i));

        % linear regression: log(err) = slope * log(h) + intercept
        p = polyfit(log_h, log_err, 1);
        slope = p(1);
        fit_vals = polyval(p, log_h);

        % plot: global truncation error vs. step size
        figure(1);
        loglog(h_avg_list(i, :), global_errors(i, :), error_data_presets{i}, ...
               'MarkerSize', 5, 'DisplayName', BT.name);
        hold on;
        loglog(10.^log_h, 10.^fit_vals, line_fit_presets{i}, ...
               'LineWidth', 2, ...
               'DisplayName', sprintf('%s Fit (slope = %.2f)', BT.name, slope));

        % second plot: global truncation error vs. number of evals
        log_num_evals = log10(num_evals_all(i, start_i:end_i));
        p_eval = polyfit(log_num_evals, log_err, 1);
        slope_eval = p_eval(1);
        fit_eval = polyval(p_eval, log_num_evals);

        figure(2);
        loglog(num_evals_all(i, :), global_errors(i, :), error_data_presets{i}, ...
               'MarkerSize', 5, 'DisplayName', sprintf('%s Data', BT.name));
        hold on;
        loglog(10.^log_num_evals, 10.^fit_eval, line_fit_presets{i}, ...
               'LineWidth', 2, ...
               'DisplayName', sprintf('%s Fit (slope = %.2f)', BT.name, slope_eval));
    end

    % plot 1: error vs. step size
    figure(1);
    title('Global Truncation Error vs. Step Size');
    xlabel('Timestep (s)');
    ylabel('Global Truncation Error');
    legend('Location', 'southeast');
    grid on;
    axis tight;
    hold off;

    % plot 2: error vs. number of function evaluations
    figure(2);
    title('Global Truncation Error vs. Number of Function Calls');
    xlabel('Number of Function Calls');
    ylabel('Global Truncation Error');
    legend('Location', 'southwest');
    grid on;
    axis tight;
    hold off;
end
