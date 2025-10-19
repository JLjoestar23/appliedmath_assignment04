% Rate function describing Newton's law of gravitation
%
% INPUTS:
% t: the current time
% V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
% orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%
% OUTPUTS:
% dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t, V, orbit_params)
    % unpack parameters
    m_s = orbit_params.m_sun;
    G = orbit_params.G;

    % unpack state variables
    x = V(1);
    y = V(2);

    % compute distance from sun
    r = sqrt(x^2 + y^2);

    % constant used in state matrix
    c = (-G*m_s)/r^3;

    % Construct derivative vector
    dVdt = [0 0 1 0; 0 0 0 1; c 0 0 0; 0 c 0 0]*V;
end
