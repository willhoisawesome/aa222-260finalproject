%% AA 260 Project
% Written by Shane Tokishi, Will Ho, Valeria Cartagena
% Last updated: 5/17/25

%% To-Do

% - think about cryogenic fuel penalty because Table 6.2 is based on statistical
% data, which does not include H2 aircraft. If we were to just add a tank
% penalty we are double counting a JetA fuel system and a H2 fuel system

%% TSFC Estimates

% Reference AA141-Propulsion Lecture 11 from Alonso
% ~ Slide 19
% Reference aircraft is the A320neo

%%
function [output, output2] = AA_260_Engine_Model(alt, Mach)
    % Operating conditions
    given.alt = alt; % cruise altitude, m
    given.alt_ft = given.alt * 3.28084; % cruise altitude, ft
    given.M = Mach; % cruise Mach number
    given.gamma = 1.4; % specific heat ratio for air
    given.R = 287; % specific gas constant for air, J/kg/K
    
    % Engine specifications
    given.TSFC0 = 0.40; % (!!!) static, sea-level installed thrust at MCT (maximum continuous thrust [normal operation])
    given.n = 0.8; % n = 0.8 is a good approximation for PW4056 (A330, B777, HBPR)
    given.ne = 2; % number of engines
    given.T0 = 26345; % PW1127G-JM max continuous sea level static thrust, lbf
    given.N1 = 0.8; % (!!!) spool speed (in percent, 0-100%) / 100
    given.m = 0.7; % m = 0.7 for h < 36151 ft; m = 1 for h >= 36151 ft
    given.a1 = -8.5e-4; % (!!!) taken from slide 20
    given.a2 = 5.5e-7; % (!!!) taken from slide 20
    given.LHV_H2 = 51500; % lower heating value of H2
    given.LHV_JetA = 18500; % lower heating value of Jet A
    
    % Determine pressure, temperature, and density ratios fror a given altitude
    atmosphere = Isentropic_Atmosphere(given.alt); % pressure, temperature, density values using isentropic atmosphere model
    a = sqrt(given.gamma * given.R * atmosphere.T); % speed of sound, m/s
    
    % Calculate TSFC and total thrust
    TSFC = zeros(size(given.M)); % initialize TSFC array
    v = zeros(size(given.M)); % initialize velocity array
    T = zeros(size(given.M)); % initialize thrust array
    
    % Rows correspond to cruise Mach, columns correspond to cruise alt
    for jj = 1:length(given.alt)
        for ii = 1:length(given.M)
            v(ii, jj) = given.M(ii) * a(jj) * 3.28084; % speed, ft/s
            TSFC(ii, jj) = given.TSFC0 * sqrt(atmosphere.t_ratio(jj)) * (1 + given.M(ii))^given.n; % thrust-specific fuel consumption, lb_fuel/hr/lbf_thrust
            T(ii, jj) = given.ne * given.N1 * given.T0 * atmosphere.r_ratio(jj)^given.m * (1 + given.a1*v(ii, jj) + given.a2*v(ii, jj)^2);
        end
    end
    
    output.T = T; % assign to output structure
    output.TSFC = TSFC; % * given.LHV_JetA / given.LHV_H2; % assign to output stucture (converted to H2 for a given thrust)
    output.mdot_f = TSFC .* T; % fuel consumption, lbm/hr
    output.v = v; % assign to output structure
    output2 = atmosphere; % assign to output stucture
    
    % figure; plot(given.M, TSFC); xlabel('M'); ylabel('TSFC [$\frac{\mathrm{lbm}}{\mathrm{hr} \cdot \mathrm{lbf}}$]'); title('TSFC vs. Cruise Mach')
    % figure; plot(given.M, T / 1000); xlabel('M'); ylabel('MCT Thrust [$\mathrm{kip}$]'); title('Maximum Continuous Thrust vs. Cruise Mach');
    % figure; plot(given.M, mdot_f); xlabel('M'); ylabel('$\dot{m}_f [\frac{\mathrm{lbm}}{\mathrm{hr}}$]'); title('Fuel Mass Flow Rate vs. Cruise Mach');
end