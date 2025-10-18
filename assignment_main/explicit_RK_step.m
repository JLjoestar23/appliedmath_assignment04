% This function computes the value of X at the next time step
% for any arbitrary RK method
% 
% INPUTS:
% rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
% t: the value of time at the current step
% XA: the value of X(t)
% h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
% BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%
% OUTPUTS:
% XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
% num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB, num_evals] = explicit_RK_step(rate_func_in, t, XA, h, BT_struct)

    % order of the chosen method
    order = length(BT_struct.C);

    % pre-allocate vector of approximations
    % length of vector should match the order
    K = zeros(order, 1);
    
    % solve for the initial K_1 approx
    K(1) = rate_func_in(t, XA);
    
    num_evals = 1; % update to 1 eval
    
    % loop through to solve each following K_i approx
    for i=2:order
        % initialize the summed a_{i,j} * k_j term
        sum_approx = 0;
        % update the summed term accordingly
        for j=1:order-1
            sum_approx = sum_approx + BT_struct.A(i, j) * K(j);
        end
        % update the current K approx
        K(i) = rate_func_in(t + c(i)*h, XA + h*sum_approx);
        num_evals = num_evals + 1; % update with additional evals
    end
    
    % initialize the summed b_i * K_i term
    lin_combo = 0;
    
    % iterate and update the summed term
    for i=1:order
        lin_combo = lin_combo + b(i)*K(i);
    end
    
    % use values to approximate the next timestep
    XB = XA + h*lin_combo;

end