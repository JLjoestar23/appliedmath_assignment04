% Runs numerical integration arbitrary RK method using variable time steps
% 
% INPUTS:
% rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
% tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
% X0: the vector describing the initial conditions, X(t_start)
% BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
% p: how error scales with step size (error = k*hˆp)
% error_desired: the desired local truncation error at each step
%
% OUTPUTS:
% t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
% X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
% h_avg: the average step size
% num_evals: total number of calls made to rate_func_in during the integration
function [t_list, X_list, h_avg, num_evals] = explicit_RK_variable_step_integration(rate_func_in, tspan, X0, BT_struct, p, error_desired)

    % Actual step size
    h_avg = (tspan(2) - tspan(1)) / N;

    % Time vector
    t_list = linspace(tspan(1), tspan(2), N+1);

    % Preallocate solution array (each row = one time step)
    X_list = zeros(length(X0), N+1);
    X_list(:, 1) = X0;

    num_evals = 0; % initialize counter
    
    for i = 1:N
        % Evaluate next step
        [XB, evals] = explicit_RK_step(rate_func_in, t_list(i), X_list(:, i), h_avg, BT_struct);
        
        % Store result as row
        X_list(:, i+1) = XB;

        % Accumulate evaluations
        num_evals = num_evals + evals;
    end
end