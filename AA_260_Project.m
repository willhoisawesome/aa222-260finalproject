%% AA 260 Project
% Written by Shane Tokishi, Will Ho, Valeria Cartagena
% Last updated: 5/8/25

%% To-Do
% - update isentropic model to return in imperial units (currently uses
% meters)
% - clean up isentropic model function
% - clean up units everywhere (mdot_f is not metric)

%% TSFC Estimates

% Reference AA141-Propulsion Lecture 11 from Alonso
% ~ Slide 19
% FILL IN WITH ASSUMPTIONS
% Reference aircraft (for now) is the A320neo

%%
clc; clear all; close all;

% Operating conditions
given.alt = 1100; % cruise altitude, m
given.M = linspace(0, 0.9); % cruise Mach number
given.gamma = 1.4; % specific heat ratio for air
given.R = 287; % specific gas constant for air, J/kg/K

% Engine specifications
given.TSFC0 = 0.40; % (!!!) static, sea-level installed thrust at MCT (maximum continuous thrust [normal operation])
given.n = 0.8; % n = 0.8 is a good approximation for PW4056 (A330, B777, HBPR)
given.ne = 2; % number of engines
given.T0 = 20305; % PW1919G max continuous sea level static thrust, lbf
given.N1 = 0.8; % (!!!) spool speed (in percent, 0-100%) / 100
given.m = 0.7; % m = 0.7 for h < 36151 ft; m = 1 for h >= 36151 ft
given.a1 = -8.5e-4; % (!!!) taken from slide 20
given.a2 = 5.5e-7; % (!!!) taken from slide 20

% Determine pressure, temperature, and density ratios fror a given altitude
atmosphere = Isentropic_Atmosphere(given.alt); % pressure, temperature, density values using isentropic atmosphere model
a = sqrt(given.gamma * given.R * atmosphere.T); % speed of sound, m/s

% Calculate TSFC and total thrust
TSFC = zeros(size(given.M)); % initialize TSFC array
v = zeros(size(given.M)); % initialize velocity array
T = zeros(size(given.M)); % initialize thrust array

for ii = 1:length(given.M)
    v(ii) = given.M(ii) * a; % speed, m/s
    TSFC(ii) = given.TSFC0 * sqrt(atmosphere.t_ratio) * (1 + given.M(ii))^given.n; % thrust-specific fuel consumption, lb_fuel/hr/lbf_thrust
    T(ii) = given.ne * given.N1 * given.T0 * atmosphere.r_ratio^given.m * (1 + given.a1*v(ii) + given.a2*v(ii)^2);
end

mdot_f = TSFC .* T; % fuel consumption, kg/s

figure; plot(given.M, TSFC); xlabel('M'); ylabel('TSFC');
figure; plot(given.M, T / 1000); xlabel('M'); ylabel('MCT Thrust, kip');
figure; plot(given.M, mdot_f); xlabel('M'); ylabel('mdot_f');