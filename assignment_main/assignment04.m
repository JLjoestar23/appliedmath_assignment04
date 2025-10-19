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

%% Local truncation error analysis

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