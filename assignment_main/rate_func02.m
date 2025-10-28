% example function 02
% describes undamped second order system
function dXdt = rate_func02(t, X)
    dXdt = [0,-1;1,0]*X;
end