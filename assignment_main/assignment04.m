%% Defining Butcher Tableau values for various RK methods

% Huen's 2nd order method Butcher Tableau
Heun2 = struct();
Heun2.C = [0, 1];
Heun2.A = [0, 0;
            1, 0];
Heun2.B = [1/2, 1/2];
Heun2.name = 'Heun2';

% VanDerHouwen/Wray's 3rd order method Butcher Tableau
VanDerHouwenWray3 = struct();
VanDerHouwenWray3.C = [0, 8/15, 2/3];
VanDerHouwenWray3.A = [0,     0,     0;
                       8/15,  0,     0;
                       1/4,   5/12,  0];
VanDerHouwenWray3.B = [1/4, 0, 3/4];
VanDerHouwenWray3.name = 'VDHWray3';

% Classic Runge-Kutta 4th order method Butcher Tableu
ClassicRK4 = struct();
ClassicRK4.C = [0, 1/2, 1/2, 1];
ClassicRK4.A = [0,   0,   0,   0;
                 1/2, 0,   0,   0;
                 0,   1/2, 0,   0;
                 0,   0,   1,   0];
ClassicRK4.B = [1/6, 1/3, 1/3, 1/6];
ClassicRK4.name = 'RK4';

%% Comparing various RK methods

% define physical/orbital parameters
orbit_params = struct();
orbit_params.m_sun = 1.9891e30; % mass of sun (kg)
orbit_params.m_planet = 5.972e24; % mass of earth (kg)
orbit_params.G = 6.6743e-11; % gravitational constant

% define rate function handle
rate_func_in = @(t,V) gravity_rate_func(t, V, orbit_params); 

% define initial conditions for a circular orbit
x0 = 1.495978707e11; % 1 astronomical unit in meters
y0 = 0;
r0 = sqrt(x0^2 + y0^2);

% calculate velocity required for circular orbit (e = 0)
% choose arbitrarily scaled down velocity for elliptical orbit
scale_factor = 0.7;
v0 = scale_factor*sqrt(orbit_params.G * orbit_params.m_sun / sqrt(x0^2 + y0^2));
vx0 = 0; % velocity perpendicular to radius
vy0 = v0; % counterclockwise motion

V0 = [x0; y0; vx0; vy0]; % initial state vector

% define time span (aim for 1 circular orbital period)
a = r0; % semi-major axis is equivalent to radius of orbit
T_orbit = 2*pi*sqrt(a^3/(orbit_params.G*orbit_params.m_sun));
tspan = [0, T_orbit];

h_ref = T_orbit/100; % want to plot 100 points to exaggerate error

% Plot orbit
method_list = {'Huen2', 'VDHWray3', 'RK4'};
BT_list = {Heun2, VanDerHouwenWray3, ClassicRK4};

figure;
hold on;
for i=1:length(BT_list)
    % call the corresponding RK method
    [t_list, V_list, ~, ~] = explicit_RK_fixed_step_integration(rate_func_in, tspan, V0, h_ref, BT_list{i});
    plot(t_list, V_list(1, :), '.-', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', method_list{i});
end
% compute analytical solution
V_analytical = compute_planetary_motion(t_list, V0, orbit_params);
plot(t_list, V_analytical(:, 1), '-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
legend();
xlabel('Time');
ylabel('X Position');
xlim([t_list(1) t_list(end)]);
title('Comparing RK Methods');
grid on;
hold off;

% plot orbit approximation
[~, V_list, ~, ~] = explicit_RK_fixed_step_integration(rate_func_in, tspan, V0, h_ref, Heun2);
figure;
hold on;
plot(V_list(1, :), V_list(2, :), '.-', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Numerical');
plot(V_analytical(:, 1), V_analytical(:, 2), '-', 'LineWidth', 1.5, 'DisplayName', 'Analytical')
plot(0, 0, 'r.', 'MarkerSize', 20, 'DisplayName', 'sun');
legend();
xlabel('X Position');
ylabel('Y Position');
title('Orbit Comparison');
grid on;
axis equal;
hold off;

%% Local truncation error analysis

% define physical/orbital parameters
orbit_params = struct();
orbit_params.m_sun = 1.9891e30; % mass of sun (kg)
orbit_params.m_planet = 5.972e24; % mass of earth (kg)
orbit_params.G = 6.6743e-11; % gravitational constant

% define initial conditions for a circular orbit
x0 = 1.495978707e11; % 1 astronomical unit in meters
y0 = 0;

% calculate velocity required for circular orbit (e = 0)
% choose arbitrarily scaled down velocity for elliptical orbit
scale_factor = 0.7;
v0 = scale_factor*sqrt(orbit_params.G * orbit_params.m_sun / sqrt(x0^2 + y0^2));
vx0 = 0; % velocity perpendicular to radius
vy0 = v0; % counterclockwise motion

V0 = [x0; y0; vx0; vy0]; % initial state vector

t_ref = 1e6; % arbitrary time
num_trials = 1000;

BT_list = {Heun2, VanDerHouwenWray3, ClassicRK4}; % define RK method list
rate_func_in = @(t,V) gravity_rate_func(t, V, orbit_params); % define rate function handle
analytical_soln = @(t) compute_planetary_motion(t, V0, orbit_params);

local_trunc_error_analysis(BT_list, rate_func_in, analytical_soln, t_ref, num_trials);

%% Global truncation error analysis

% define physical/orbital parameters
orbit_params = struct();
orbit_params.m_sun = 1.9891e30; % mass of sun (kg)
orbit_params.m_planet = 5.972e24; % mass of earth (kg)
orbit_params.G = 6.6743e-11; % gravitational constant

% define initial conditions for a circular orbit
x0 = 1.495978707e11; % 1 astronomical unit in meters
y0 = 0;
r0 = sqrt(x0^2 + y0^2);

% calculate velocity required for circular orbit (e = 0)
% choose arbitrarily scaled down velocity for elliptical orbit
scale_factor = 0.7;
v0 = scale_factor*sqrt(orbit_params.G * orbit_params.m_sun / sqrt(x0^2 + y0^2));
vx0 = 0; % velocity perpendicular to radius
vy0 = v0; % counterclockwise motion

V0 = [x0; y0; vx0; vy0]; % initial state vector

% define time span (semi-arbirtarily)
a = r0; % semi-major axis is equivalent to radius of orbit
T_orbit = pi*sqrt(a^3/(orbit_params.G*orbit_params.m_sun));

tspan = [0, T_orbit];

BT_list = {Heun2, VanDerHouwenWray3, ClassicRK4}; % define RK method list
rate_func_in = @(t,V) gravity_rate_func(t, V, orbit_params); % define rate function handle
analytical_soln = @(t) compute_planetary_motion(t, V0, orbit_params);

global_trunc_error_analysis(BT_list, rate_func_in, analytical_soln, tspan, 250);

%% Conservation of physical quantities

% define physical/orbital parameters
orbit_params = struct();
orbit_params.m_sun = 1.9891e30; % mass of sun (kg)
orbit_params.m_planet = 5.972e24; % mass of earth (kg)
orbit_params.G = 6.6743e-11; % gravitational constant

% put them into variables for ease of use
m_sun = orbit_params.m_sun;
m_planet = orbit_params.m_planet;
G = orbit_params.G;

% define function handles
rate_func_in = @(t,V) gravity_rate_func(t, V, orbit_params);
step_func = @leapfrog_step;

% define initial conditions for a circular orbit
x0 = 1.495978707e11; % 1 astronomical unit in meters
y0 = 0;
r0 = sqrt(x0^2 + y0^2);

% calculate velocity required for circular orbit (e = 0)
% choose arbitrarily scaled down velocity for elliptical orbit
scale_factor = 0.7;
v0 = scale_factor*sqrt(orbit_params.G * orbit_params.m_sun / sqrt(x0^2 + y0^2));
vx0 = 0; % velocity perpendicular to radius
vy0 = v0; % counterclockwise motion

V0 = [x0; y0; vx0; vy0]; % initial state vector

% define time span (aim for 1.5 circular orbital period)
a = r0; % semi-major axis is equivalent to radius of orbit
T_orbit = 3*pi*sqrt(a^3/(orbit_params.G*orbit_params.m_sun));
tspan = [0, T_orbit];

h_ref = T_orbit/1000;

% Plot orbit
method_list = {'2nd Order', '3rd Order', '4th Order'};
BT_list = {Heun2, VanDerHouwenWray3, ClassicRK4};

% plot comparing conservation of mechanical energy
figure();
for i=1:length(BT_list)
    % call the corresponding RK method
    [t_list, V_list, ~, ~] = explicit_RK_fixed_step_integration(rate_func_in, tspan, V0, h_ref, BT_list{i});
    
    x = V_list(1, :);
    y = V_list(2, :);
    x_dot = V_list(3, :);
    y_dot = V_list(4, :);
    
    % calculate approximated E over time
    E = 0.5 * m_planet * (x_dot.^2 + y_dot.^2) - (m_sun * m_planet * G)./sqrt(x.^2 + y.^2);
    % initial state E
    E_i = 0.5 * m_planet * (x_dot(1).^2 + y_dot(1).^2) - (m_sun * m_planet * G)./sqrt(x(1).^2 + y(1).^2);
    % normalized error in E
    E_err = (E - E_i) / abs(E_i);
    plot(t_list, E_err, '.', 'MarkerSize', 8, 'DisplayName', method_list{i});
    hold on;
end

% call leapfrog method
[t_list, V_list, ~, ~] = fixed_step_integration(rate_func_in, step_func, tspan, V0, h_ref);
V_list = V_list';
x = V_list(1, :);
y = V_list(2, :);
x_dot = V_list(3, :);
y_dot = V_list(4, :);
% calculate approximated E over time
E = 0.5 * m_planet * (x_dot.^2 + y_dot.^2) - (m_sun * m_planet * G)./sqrt(x.^2 + y.^2);
% initial state E
E_i = 0.5 * m_planet * (x_dot(1).^2 + y_dot(1).^2) - (m_sun * m_planet * G)./sqrt(x(1).^2 + y(1).^2);
% normalized error in E
E_err = (E - E_i) / abs(E_i);
plot(t_list, E_err, '.', 'MarkerSize', 8, 'DisplayName', 'Leapfrog');

% plot aesthetics
legend('Location', 'northwest');
xlabel('Time');
ylabel('Normalized Mechanical Energy Error');
xlim([t_list(1) t_list(end)]);
title('Conservation of E between RK Methods');
grid on;
hold off;

% plot comparing conservation of angular momentum
figure();
for i=1:length(BT_list)
    % call the corresponding RK method
    [t_list, V_list, ~, ~] = explicit_RK_fixed_step_integration(rate_func_in, tspan, V0, h_ref, BT_list{i});
    
    x = V_list(1, :);
    y = V_list(2, :);
    x_dot = V_list(3, :);
    y_dot = V_list(4, :);
    
    % calculate approximated H over time
    H = m_planet * (x .* y_dot - y .* x_dot);
    % initial state E
    H_i = m_planet * (x(1) .* y_dot(1) - y(1) .* x_dot(1));
    % normalized error in E
    H_err = (H - H_i) / abs(H_i);
    plot(t_list, H_err, '.', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', method_list{i});
    hold on;
end

% call leapfrog method
[t_list, V_list, ~, ~] = fixed_step_integration(rate_func_in, step_func, tspan, V0, h_ref);
V_list = V_list';
x = V_list(1, :);
y = V_list(2, :);
x_dot = V_list(3, :);
y_dot = V_list(4, :);
% calculate approximated H over time
H = m_planet * (x .* y_dot - y .* x_dot);
% initial state E
H_i = m_planet * (x(1) .* y_dot(1) - y(1) .* x_dot(1));
% normalized error in E
H_err = (H - H_i) / abs(H_i);
plot(t_list, H_err, '.', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Leapfrog');

% plot aesthetics
legend('Location', 'northwest');
xlabel('Time');
ylabel('Normalized Angular Momentum Error');
xlim([t_list(1) t_list(end)]);
title('Conservation of H between RK Methods');
grid on;
hold off;

%% Analysis on differently ordered next step estimates

% Dormand Prince Butcher Tableau
DormandPrince = struct();
DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
                    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
DormandPrince.A = [0,0,0,0,0,0,0;
                    1/5, 0, 0, 0,0,0,0;...
                    3/40, 9/40, 0, 0, 0, 0,0;...
                    44/45, -56/15, 32/9, 0, 0, 0,0;...
                    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
                    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
                    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

% initializing function handles
rate_func_in = @rate_func02;
analytical_soln = @solution02;

% number of points to plot
num_trials = 1000;

% arbitrary time
t = 0.314;

% initialize step sizes (log-spaced)
h_list = logspace(-6, 1, num_trials); % seconds

% initialize storage
loc_trunc_err_1 = zeros(1, num_trials);
loc_trunc_err_2 = zeros(1, num_trials);
order_diff = zeros(1, num_trials);
local_h_diff = zeros(1, num_trials);

% plot style presets
error_data_presets = {'r.', 'b.', 'm-', 'k-'};

% compute local truncation errors for this method
for i = 1:num_trials
    % get current state from the analytical solution
    XA = analytical_soln(t);

    % get next state from the analytical solution
    XB_ana = analytical_soln(t + h_list(i));

    % compute both embedded next timestep estimates
    [XB1, XB2, ~] = embedded_RK_step(rate_func_in, t, XA, h_list(i), DormandPrince);
    
    % compute both local truncation errors
    loc_trunc_err_1(i) = norm(XB1 - XB_ana);
    loc_trunc_err_2(i) = norm(XB2 - XB_ana);

    % compute abs diff between 2 estimates
    order_diff(i) = norm(loc_trunc_err_1 - loc_trunc_err_2);

    % Compute difference in the analytical trajectory over one step
    local_h_diff(i) = norm(XB_ana - XA);
end

% plot data as a function of h
figure();
loglog(h_list, loc_trunc_err_1, 'r.', 'MarkerSize', 5, 'DisplayName', 'XB1 error');
hold on;
loglog(h_list, loc_trunc_err_2, 'b.', 'MarkerSize', 5, 'DisplayName', 'XB2 error');
loglog(h_list, order_diff, 'm--', 'LineWidth', 2, 'DisplayName', '|XB1-XB2|');
loglog(h_list, local_h_diff, 'k--', 'LineWidth', 2, 'DisplayName', 'timestep difference');

% labels and aesthetics
title('Local Truncation Error vs. Step Size');
xlabel('Timestep (s)');
ylabel('Local Truncation Error');
legend('Location', 'northwest');
grid on;
axis tight;
hold off;

% finding slopes of loc trunc err
start= 0.6 * num_trials;
fin = 0.8 * num_trials;
[p1, ~] = polyfit(log10(h_list(start:fin)), log10(loc_trunc_err_1(start:fin)), 1);
[p2, ~] = polyfit(log10(h_list(start:fin)), log10(loc_trunc_err_2(start:fin)), 1);
fprintf('Estimated Local Truncation Error Order: XB1 = %.2f, XB2 = %.2f\n', p1(1), p2(1));

% plot data as a function of |XB1-XB2|
figure();
loglog(order_diff, loc_trunc_err_1, 'r.', 'MarkerSize', 5, 'DisplayName', 'XB1 error');
hold on;
loglog(order_diff, loc_trunc_err_2, 'b.', 'MarkerSize', 5, 'DisplayName', 'XB2 error');

% labels and aesthetics
title('Local Truncation Error vs. Order Difference');
xlabel('|XB1-XB2|');
ylabel('Local Truncation Error');
legend('Location', 'northwest');
grid on;
axis square;
hold off;
hold off;

%% Deliverables

% Dormand Prince Butcher Tableau
DormandPrince = struct();
DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
                    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
DormandPrince.A = [0,0,0,0,0,0,0;
                    1/5, 0, 0, 0,0,0,0;...
                    3/40, 9/40, 0, 0, 0, 0,0;...
                    44/45, -56/15, 32/9, 0, 0, 0,0;...
                    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
                    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
                    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

% define physical/orbital parameters
orbit_params = struct();
orbit_params.m_sun = 1.9891e30; % mass of sun (kg)
orbit_params.m_planet = 5.972e24; % mass of earth (kg)
orbit_params.G = 6.6743e-11; % gravitational constant
m_sun = orbit_params.m_sun;
m_planet = orbit_params.m_planet;
G = orbit_params.G;

% define rate function handle
rate_func_in = @(t,V) gravity_rate_func(t, V, orbit_params);
analytical_soln = @(t) compute_planetary_motion(t, V0, orbit_params);

% define initial conditions for a circular orbit
x0 = 1.495978707e11; % 1 astronomical unit in meters
y0 = 0;
r0 = sqrt(x0^2 + y0^2);

% calculate velocity required for circular orbit (e = 0)
% choose arbitrarily scaled down velocity for elliptical orbit
scale_factor = 0.7;
v0 = scale_factor*sqrt(orbit_params.G * orbit_params.m_sun / sqrt(x0^2 + y0^2));
vx0 = 0; % velocity perpendicular to radius
vy0 = v0; % counterclockwise motion

V0 = [x0; y0; vx0; vy0]; % initial state vector

% define time span (aim for 1 circular orbital period)
a = r0; % semi-major axis is equivalent to radius of orbit
T_orbit = pi*sqrt(a^3/(orbit_params.G*orbit_params.m_sun));
tspan = [0, T_orbit];

% number of timesteps in parameter sweep
num_trials = 250;

% generate log-spaced step sizes
desired_error_list = logspace(-8, 8, num_trials);
desired_h_list = logspace(3, 7, num_trials);

% initialize storage arrays
global_trunc_errors = zeros(2, num_trials);
h_avg_list = zeros(2, num_trials);
num_evals = zeros(2, num_trials);
failed_fracs = zeros(1, num_trials);

t_ref = linspace(tspan(1), tspan(2), 1e5);
V_ref = compute_planetary_motion(t_ref, V0, orbit_params)';
X_final = V_ref(:, end);
X0 = V_ref(:, 1);

% compute global truncation error for each step size
for i = 1:num_trials

    % Copmute numerical solution using variable step method
    [~, X_list_num_var, h_avg_var, num_evals_var, failed_frac_var] = ...
    explicit_RK_variable_step_integration(rate_func_in, tspan, X0, DormandPrince, 5, desired_error_list(i));

    [~, X_list_num_fix, h_avg_fix, num_evals_fix] = ...
    explicit_RK_fixed_step_integration(rate_func_in, tspan, X0, desired_h_list(i), DormandPrince);

    % Compute relevant analysis values
    global_trunc_errors(1, i) = norm(X_list_num_var(:, end) - X_final);
    h_avg_list(1, i) = h_avg_var;
    num_evals(1, i) = num_evals_var;
    global_trunc_errors(2, i) = norm(X_list_num_fix(:, end) - X_final);
    h_avg_list(2, i) = h_avg_fix;
    num_evals(2, i) = num_evals_fix;
    failed_fracs(i) = failed_frac_var;
end

%% Plots for above

figure();
loglog(h_avg_list(1, :), global_trunc_errors(1, :), 'r.', 'MarkerSize', 8, 'DisplayName', 'Variable h');
hold on;
loglog(h_avg_list(2, :), global_trunc_errors(2, :), 'b.', 'MarkerSize', 8, 'DisplayName', 'Fixed h');
legend();
title('Global Truncation Error vs. Average Step Size');
xlabel('Timestep (s)');
ylabel('Global Truncation Error');
legend('Location', 'northwest');
grid on;
axis tight;
hold off;

figure();
loglog(num_evals(1, :), global_trunc_errors(1, :), 'r.', 'MarkerSize', 8, 'DisplayName', 'Variable h');
hold on;
loglog(num_evals(2, :), global_trunc_errors(2, :), 'b.', 'MarkerSize', 8, 'DisplayName', 'Fixed h');
legend();
title('Global Truncation Error vs. # of Evaluations');
xlabel('# of Evaluations');
ylabel('Global Truncation Error');
legend('Location', 'northeast');
grid on;
axis tight;
hold off;

figure();
semilogx(h_avg_list(1, :), failed_fracs, 'r.', 'MarkerSize', 8, 'DisplayName', 'Failure Rate');
hold on;
legend();
title('Failure Rate vs Average Step Size');
xlabel('Timestep (s)');
ylabel('Failure Rate');
legend('Location', 'northeast');
grid on;
axis tight;
hold off;

%% Plotting orbit

% Dormand Prince Butcher Tableau
DormandPrince = struct();
DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
                    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
DormandPrince.A = [0,0,0,0,0,0,0;
                    1/5, 0, 0, 0,0,0,0;...
                    3/40, 9/40, 0, 0, 0, 0,0;...
                    44/45, -56/15, 32/9, 0, 0, 0,0;...
                    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
                    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
                    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

% define physical/orbital parameters
orbit_params = struct();
orbit_params.m_sun = 1.9891e30; % mass of sun (kg)
orbit_params.m_planet = 5.972e24; % mass of earth (kg)
orbit_params.G = 6.6743e-11; % gravitational constant
m_sun = orbit_params.m_sun;
m_planet = orbit_params.m_planet;
G = orbit_params.G;

% define rate function handle
rate_func_in = @(t,V) gravity_rate_func(t, V, orbit_params); 

% define initial conditions for a circular orbit
x0 = 1.495978707e11; % 1 astronomical unit in meters
y0 = 0;
r0 = sqrt(x0^2 + y0^2);

% calculate velocity required for circular orbit (e = 0)
% choose arbitrarily scaled down velocity for elliptical orbit
scale_factor = 0.3;
v0 = scale_factor*sqrt(orbit_params.G * orbit_params.m_sun / sqrt(x0^2 + y0^2));
vx0 = 0; % velocity perpendicular to radius
vy0 = v0; % counterclockwise motion

V0 = [x0; y0; vx0; vy0]; % initial state vector

% define time span (aim for 1 circular orbital period)
a = r0; % semi-major axis is equivalent to radius of orbit
T_orbit = 0.75*pi*sqrt(a^3/(orbit_params.G*orbit_params.m_sun));
tspan = [0, T_orbit];

% compute analytical soln
V_analytical = compute_planetary_motion(linspace(tspan(1), tspan(2), 1000), V0, orbit_params);

% call the corresponding RK method
[t_list, V_list, h_avg, ~, ~] = explicit_RK_variable_step_integration(rate_func_in, tspan, V0, DormandPrince, 5, 1e2);

% plot orbit approximation
figure;
hold on;
plot(V_list(1, :), V_list(2, :), 'ro-', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', 2, 'DisplayName', 'planet');
plot(0, 0, 'r.', 'MarkerSize', 40, 'DisplayName', 'sun');
legend();
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Orbit Plot');
grid on;
axis equal
hold off;

%% Step size as a function of orbit radius

r = sqrt(V_list(1, :).^2 + V_list(2, :).^2);
h_list = diff(t_list);

figure();
loglog(h_list, r(1:end-1), '.', 'MarkerSize', 5)
xlabel('Timestep Size (s)');
ylabel('Orbit Radius (m)');
title('Timestep Size vs Orbit Radius');
grid on;
axis tight;
hold off;