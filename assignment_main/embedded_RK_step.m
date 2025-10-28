% This function computes the value of X at the next time step
% for any arbitrary embedded RK method
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
% XB1: the approximate value for X(t+h) using the first row of the Tableau
% XB2: the approximate value for X(t+h) using the second row of the Tableau
% num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB1, XB2, num_evals] = embedded_RK_step(rate_func_in, t, XA, h, BT_struct)
    % parse the struct
    A = BT_struct.A;
    B = BT_struct.B;
    C = BT_struct.C;
    
    % # of stages = order of the method
    stages = length(C);

    % state dimension
    n = length(XA);

    % pre-allocate vector of approximations
    % length of vector should match the order
    K = zeros(n, stages);
    
    % compute all K_i stages
    for i=1:stages
        % evaluate the sum of a_{i,j}*k_j terms as dot product
        sum_val = K*(A(i, :)');
        % evaluate the ith K approx
        K(:, i) = rate_func_in(t + C(i)*h, XA + h*(sum_val));
    end
    
    % compute both embedded next timestep estimates
    XB1 = XA + h * (K*B(1,:)');
    XB2 = XA + h * (K*B(2,:)');

    % # of evals should equal the # number of stages
    num_evals = stages;
end