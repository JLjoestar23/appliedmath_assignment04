function [XB, num_evals] = leapfrog_step(rate_func_in, t, XA, h)
    
    r = XA(1:2);
    v = XA(3:4);

    current_rate = rate_func_in(t, XA);
    a = current_rate(3:4);

    v_half = v + 0.5 * a * h;

    r_next = r + v_half * h;
    
    next_temp_state = [r_next, v_half];
    next_rate = rate_func_in(t + h, next_temp_state);
    v_next = v_half + 0.5 * next_rate * h;

    XB = [r_next; v_next];

    num_evals = 2;
end