function [XB, num_evals] = leapfrog_step(rate_func_in, t, XA, h)
    % Extract position and velocity
    r = XA(1:2);
    v = XA(3:4);

    % Compute acceleration at current state
    current_rate = rate_func_in(t, XA);
    a = current_rate(3:4);

    % Half-step velocity
    v_half = v + 0.5 * a * h;

    % Full-step position
    r_next = r + v_half * h;

    % Compute acceleration at next position
    next_temp_state = [r_next; v_half];
    next_rate = rate_func_in(t + h, next_temp_state);
    a_next = next_rate(3:4);

    % Complete velocity step
    v_next = v_half + 0.5 * a_next * h;

    % Combine into next state
    XB = [r_next; v_next];

    % Count of rate function evaluations
    num_evals = 2;
end
