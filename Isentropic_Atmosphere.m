%% Atmospheric conditions as a function of altitude
% Written by Shane Tokishi
% Last updated: 5/9/25

% Assumptions: 
% Code assumes isentropic, adiabatic, and significant buoyancy forces:
% (air density is large enough that resulting buoyancy forces due to 
% temperature differences are large enough to move infinitesimally small 
% parcels of air to its equilibrium altitude)

% These assumptions are valid throughout the troposphere and into the
% stratosphere, up to an altitude of about 18 km according to Dahm

% Past 18 km, the density becomes too low and thus buoyancy forces are too
% small, and solar UV absorbtion causes the adiabatic assumption to be poor

function output = Isentropic_Atmosphere(alt)
%% Constants
given.gamma = 1.4; % air
given.p_surf = 101325; % N/m^2
given.rho_surf = 1.23; % kg/m^3
given.g = 9.81; % m/s
given.T_surf = 288; % K
given.R = 287; % J/kg/K
given.z_star = given.R * given.T_surf / given.g; % m

%% Pressure
p_plot.p_fun = @(z) ((1-(((given.gamma-1)./given.gamma)...
    .* (z./given.z_star))).^(given.gamma ./ (given.gamma-1))) * given.p_surf; % local pressure function, Pa

output.p = p_plot.p_fun(alt); % local pressure, Pa
output.p_ratio = p_plot.p_fun(alt) / given.p_surf; % pressure ratio to sea level

%% Temperature
t_plot.t_fun = @(z) (1-(((given.gamma-1)./given.gamma)...
    .* (z./given.z_star))) * given.T_surf; % local temperature function, K

output.T = t_plot.t_fun(alt); % local temperature, K
output.t_ratio = output.T / given.T_surf; % temperature ratio to sea level

%% Density
r_fun = @(z) ((1-(((given.gamma-1)./given.gamma)...
    .* (z./given.z_star))).^(1 ./ (given.gamma-1))) * given.rho_surf; % local density function, kg/m^3

output.rho = r_fun(alt); % local density, kg/m^3
output.r_ratio = output.rho / given.rho_surf; % density ratio to sea level
end