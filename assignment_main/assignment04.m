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
method_list = {'2nd Order', '3rd Order', '4th Order'};
BT_list = {Heun2, VanDerHouwenWray3, ClassicRK4};

figure;
hold on;
for i=1:length(BT_list)
    % call the corresponding RK method
    [t_list, V_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, tspan, V0, h_ref, BT_list{i});
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
[t_list, V_list, ~, ~] = explicit_RK_fixed_step_integration(rate_func_in, tspan, V0, h_ref, Heun2);
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

global_trunc_error_analysis(BT_list, rate_func_in, analytical_soln, tspan, num_trials);

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
    [t_list, V_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, tspan, V0, h_ref, BT_list{i});
    
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
    plot(t_list, E_err, '.', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', method_list{i});
    hold on;
end
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
    [t_list, V_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, tspan, V0, h_ref, BT_list{i});
    
    x = V_list(1, :);
    y = V_list(2, :);
    x_dot = V_list(3, :);
    y_dot = V_list(4, :);
    
    % calculate approximated E over time
    H = m_planet * (x .* y_dot - y .* x_dot);
    % initial state E
    H_i = m_planet * (x(1) .* y_dot(1) - y(1) .* x_dot(1));
    % normalized error in E
    H_err = (H - H_i) / abs(H_i);
    plot(t_list, H_err, '.', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', method_list{i});
    hold on;
end
legend('Location', 'northwest');
xlabel('Time');
ylabel('Normalized Angular Momentum Error');
xlim([t_list(1) t_list(end)]);
title('Conservation of H between RK Methods');
grid on;
hold off;

%% Test leap frog integration

% define physical/orbital parameters
orbit_params = struct();
orbit_params.m_sun = 1.9891e30; % mass of sun (kg)
orbit_params.m_planet = 5.972e24; % mass of earth (kg)
orbit_params.G = 6.6743e-11; % gravitational constant

% put them into variables for ease of use
m_sun = orbit_params.m_sun;
m_planet = orbit_params.m_planet;
G = orbit_params.G;

% define rate function handle
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

[t_list, V_list, h_avg, num_evals] = fixed_step_integration(rate_func_in, step_func, tspan, V0, h_ref);
plot(t_list, V_list(1, :), '.-', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Leapfrog');