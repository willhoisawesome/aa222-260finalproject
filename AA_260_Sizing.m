%% AA 260 Sizing
% Written by Shane Tokishi

% This MATLAB code uses the parametric statistical approach outlined in
% table 6.1 of Daniel P. Raymer's Aircraft Design: A Conceptual Appraoch 
% (6th Edition) in order to roughly estimate initial sizing for a hydrogen
% aircraft baselined on the A320neo for AA 260 Sustainable Aviation.

% Units are variable and are labelled in the comments
% Final weights are in lbs

% -------------------------------------------------------------------------

%% Assumptions
% Parameters based on a 'General Aviation - twin engine' in Table 6.1
clc; clear; close all;

given.a = -0.90;
given.b = 1.36;
given.C1 = -0.10;
given.C2 = 0.08;
given.C3 = 0.05;
given.C4 = -0.05;
given.C5 = 0.20;

given.lambda = 27; % A320neo leading edge sweep angle
given.AR = 10.5; % A320neo aspect ratio

given.M_cruise = 0.6:0.1:0.9; % A320neo cruise mach range
given.hp_W0 = 0.036 * given.v_max^0.32; % power loading for General Aviation - twin engine (Table 5.4, pg. 120)

given.alt = 9000:1000:12000; % cruise altitude m
given.range = 2000; % nautical miles
given.W_S = [50:10:200]; % estimated wing loading, lbs/ft^2 (77 kg/m^2 is the F22 wing W/S)

%% Generate engine data over the altitude and Mach range
[turbofan, atmosphere] = AA_260_Engine_Model(given.alt, given.M_cruise);

%% Aircraft Aerodynamics (!!!) NEEDS TO BE UPDATED
% Wing loading currently done in a loop, cruise Mach and cruise alt are
% matrix operations, not sure which is better. Can probably change indices
% to be (ii, :, :) to keep iterations over Mach and alt to be matrix
% multiplication as well as wing loading as a loop

for ii = 1:length(given.W_S)

    given.coef1 = (given.b * given.AR^given.C2 * given.hp_W0^given.C3 * given.W_S(ii)^given.C4 * given.v_max^given.C5); % used in empty weight equation from Table 6.2 (pg. 149)
    given.coef2 = given.a; % used in empty weight equation from Table 6.2 (pg. 149)
    drag.cd0 = 0.022; % A320neo typical cd0
    ii = 9; % just need this so that it runs
    
    %% Fuel Weight Fractions
    % 1     Takeoff
    % 2     Climb
    % 3     Cruise
    % 10    Descent
    % 11    Landing

    % Takeoff
    weight_frac.w01 = 0.98; % taken from pg. 150 (W1/W0 ~ 0.97-0.99)

    % Climb (M0.1 --> M_cruise)
    weight_frac.w12 = 0.991 - 0.007*given.M_cruise - 0.01*given.M_cruise.^2;

    % Cruise (M_cruise)
    cruise.R = given.range * 6076; % range (nautical miles --> ft)
    cruise.C = turbofan.TSFC / 3600; % TSFC from AA_260_Engine_Model.m in lbm/hr/lbf converted to lbm/s/lbf;
    cruise.V = turbofan.v; % cruise velocities in ft/s
    cruise.q = (1/2) * atmosphere.rho / 515.4 .* cruise.V.^2; % dynamic pressure (kg/m^3 --> slug/ft^3 from Isentropic_Atmosphere.m)
    cruise.e = 0.82;
    cruise.L_D = 1 ./ ((cruise.q * drag.cd0 ./ given.W_S(ii)) + (given.W_S(ii) ./ (cruise.q * pi * given.AR * cruise.e))); % eqn. 6.13 (pg. 151)
    weight_frac.w23 = exp(-cruise.R*cruise.C / (cruise.V*cruise.L_D));
    
    % Descent
    weight_frac.w910 = 0.992; % eqn. 6.22 (pg. 153)

    % Landing
    weight_frac.w1011 = 0.995; % eqn. 6.23 (pg. 153)

    weight_frac.final = weight_frac.w01 * weight_frac.w12 * weight_frac.w23...
        * weight_frac.w910 * weight_frac.w1011; % dry mass / wet mass

    %% Iterative Sizing Process
    clear sizing

    sizing.wf_w0 = 1.06 * (1- weight_frac.final); % fuel fraction (6% accounts for reserve fuel and stuff)
    sizing.w_crew_payload = 220; % lbs (!!!) NEEDS TO BE UPDATED
    sizing.diff = 1;
    sizing.count = 1;

    while sizing.diff > 1E-3

        if sizing.count == 1
            sizing.w0_guess(sizing.count) = 1200;                
        else
            sizing.w0_guess(sizing.count) = 0.2 * (sizing.w0_guess(sizing.count - 1) - sizing.w0_calc(sizing.count-1)) + sizing.w0_calc(sizing.count-1); % ref. example on pg. 916
        end
        
        sizing.w_empty(sizing.count) = ((given.coef1 * sizing.w0_guess(1)^(-0.1) + given.coef2) * sizing.w0_guess(1)); % ref. example on pg. 875
        sizing.w_fuel(sizing.count) = sizing.wf_w0 * sizing.w0_guess(sizing.count);
        sizing.w0_calc(sizing.count) = sizing.w_fuel(sizing.count) + sizing.w_crew_payload + sizing.w_empty(sizing.count);
        sizing.diff = abs((sizing.w0_guess(sizing.count) - sizing.w0_calc(sizing.count)) / sizing.w0_guess(sizing.count));
        sizing.count = sizing.count + 1;

    end

    w0(ii) = sizing.w0_calc(end);
end

%% Total Drag (!!!) NEEDS TO BE UPDATED
sizing.s = w0 ./ given.W_S; % wing area, ft^2 - calculated based on MINIMUM wing loading
sizing.b = sqrt(given.AR .* sizing.s); % calculate wing span, ft
drag.total_lbf = drag.cd0 .* sizing.s .* cruise.q + (w0 ./ (cruise.q .* pi .* sizing.b.^2 .* cruise.e)); % total drag, lbf
drag.total_N = drag.total_lbf * 4.44822; % total drag, N
drag.T_W_max = drag.total_lbf ./ w0; % maximum thrust-to-weight ratio
drag.T_D = ACE_cruise.turbofan(6) ./ drag.total_N; % cruise thrust divided by cruise drag (should be greater than 1)

%% Trim Data (!!!) NEEDS TO BE UPDATED
sizing.w0_trim = w0(w0<150000); % trim extreme initial weight data
sizing.w0_trim_index = find(w0<150000); % find locations of valid data
given.W_S_trim = given.W_S(sizing.w0_trim_index(1):sizing.w0_trim_index(end)); % trim wing loadings
drag.T_D_trim = drag.T_D(sizing.w0_trim_index(1):sizing.w0_trim_index(end)); % trim thrust-to-drag
sizing.s_trim = sizing.s(sizing.w0_trim_index(1):sizing.w0_trim_index(end)); % trim planform area
sizing.b_trim = sizing.b(sizing.w0_trim_index(1):sizing.w0_trim_index(end)); % trim wing span

%% Plot
ThesisPlot(1, jj, given.W_S_trim, sizing.w0_trim/1000, "Wing Loading, lbs/ft^2", "Initial Weight, kip", given.alt(jj), drag.T_D_trim);
ThesisPlot(2, jj, given.W_S_trim, sizing.b_trim, "Wing Loading, lbs/ft^2", "Wing Span, ft", given.alt(jj), drag.T_D_trim);
ThesisPlot(3, jj, given.W_S_trim, sizing.s_trim, "Wing Loading, lbs/ft^2", "Wing Area, ft^2", given.alt(jj), drag.T_D_trim);    

%% Record Plotted Data
field = sprintf('alt%i', given.alt(jj)/1000);

test(jj).w0 = w0;
test(jj).s = sizing.s;
test(jj).b = sizing.b;

figure(1); sgtitle("Initial Weight vs. Wing Loading");
figure(2); sgtitle("Wing Span vs. Wing Loading");
figure(3); sgtitle("Wing Area vs. Wing Loading");