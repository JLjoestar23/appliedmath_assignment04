% Calculates the local truncation error for a generic RK method
%
% INPUTS:
% rate_func_in:  handle to rate function, f(t, X)
% analytical_soln: handle to the exact analytical solution X(t)
% t:  current time
% h:  step size
% BT_struct: struct with fields A, B, C (Butcher tableau)
%
% OUTPUTS:
% local_trunc_error: norm of (numerical - analytical) local error
% local_h_diff: norm of (exact step size change)
function [local_trunc_error, local_h_diff] = calc_local_trunc_error(BT_struct, rate_func_in, analytical_soln, t, h)
    % get current state from the analytical solution
    XA = analytical_soln(t);
    
    % get next state from the analytical solution
    XB_ana = analytical_soln(t + h);

    % compute next step approx using given RK method
    [XB_num, ~] = explicit_RK_step(rate_func_in, t, XA, h, BT_struct);

    % compute local truncation error
    local_trunc_error = norm(XB_num - XB_ana);

    % Compute difference in the analytical trajectory over one step
    local_h_diff = norm(XB_ana - XA);
end
