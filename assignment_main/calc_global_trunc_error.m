function [global_trunc_error, h_avg, num_evals] = calc_global_trunc_error(solver, rate_func_in, analytical_soln, tspan, h_ref)

    % Get the current state from the analytical solution
    X0 = analytical_soln(tspan(1)); % start from exact initial condition

    % case switching to handle multiple solver types
    switch solver
        case 'ForwardEuler'
           step_func = @forward_euler_step;

        case 'ExplicitMidpoint'
           step_func = @explicit_midpoint_step;

        case 'BackwardEuler'
           step_func = @backward_euler_step;

        case 'ImplicitMidpoint'
           step_func = @implicit_midpoint_step;

        otherwise
            warning('Invalid method');
            global_trunc_error = NaN;
            return;
    end

    % Compute next numerical step
    [t_list, X_list_num, h_avg, num_evals] = fixed_step_integration(rate_func_in, step_func, tspan, X0, h_ref);
    X_list_ana = analytical_soln(t_list(end));
    global_trunc_error = norm(X_list_num(end, :)' - X_list_ana);
end
