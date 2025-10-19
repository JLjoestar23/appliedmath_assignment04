function [global_trunc_error, h_avg, num_evals] = calc_global_trunc_error(BT_struct, rate_func_in, analytical_soln, tspan, h_ref)
    
    % Get the current state from the analytical solution
    X0 = analytical_soln(tspan(1)); % start from exact initial condition
    
    % Compute next numerical step
    [t_list, X_list_num, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, tspan, X0, h_ref, BT_struct);
    X_list_ana = analytical_soln(t_list(end));
    global_trunc_error = norm(X_list_num(:, end) - X_list_ana);
end
