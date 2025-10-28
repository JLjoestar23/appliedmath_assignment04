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
    
    % parse the struct
    A = BT_struct.A;
    B = BT_struct.B;
    C = BT_struct.C;
    
    % order of the chosen method
    stages = length(BT_struct.C);
    n = length(XA);

    % pre-allocate vector of approximations
    % length of vector should match the order
    K = zeros(n, stages);
    
    % compute all K_i stages
    for i=1:stages
        % evaluate the sum of a_{i,j}*k_j terms as dot product
        sum_val1 = K*(A(i, :)');
        % evaluate the ith K approx
        K(:, i) = rate_func_in(t + C(i)*h, XA + h*(sum_val1));
    end

    % evaluate the sum of b_i*k_i terms
    sum_val2 = K*B';

    % evaluate the next timestep estimate
    XB = XA + h*sum_val2;

    % number of evals should equal the number of stages
    num_evals = stages;

end