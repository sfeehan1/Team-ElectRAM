% Reference Solution - Homework 5 - Critical Mach Numbers
% For Aerodynamics teams
%% --------------------------------------------------
% Aircraft Parameters
alpha = linspace(-15, 15, 100); % AOA array (deg)
alpha = alpha';
AR = 3; % Aspect ratio
Swet = 1371; % Wetted area [ft^2]
S = 300; % Wing planform area [ft^2]
Lambda_LE = 40; % Leading edge sweep [deg]
NACA = 1404; % 4-digit airfoil
Length = 48.30; % Fuselage length [ft]
Amax = 25.11; % Max fuselage cross-sectional area [ft^2]
Cfe = 0.0037; % Skin friction drag coefficient (from Brandt's work as well as
% Roskam, Air Vehicle Design, Volume 1, Table 3.5
e_osw = 0.914; % Oswald span efficiency factor (dimensionless)
E_wd = 2.2; % I have no clue what this is lmfao (Found in Brandt's work)

% Additional parameters
h_ft = 35000; % Altitude in feet
M = 0.50; % Mach number
% W_lb = 31377; % Aircraft weight in lbf
Wto_S = 104.59;
W_Wto = 0.899666963;
W_S = Wto_S * W_Wto; % Aircraft wing loading in lbf/ft^2

% Obtaining t/c
tc = 0.04

% Derived Mach points
M_CD0_max = real(1 / (cosd(Lambda_LE)^0.2)); % Determines the Mach at CD0 max or something?
M_Crit = 1.0 - 0.065 * (cosd(Lambda_LE) * 100 * (tc))^0.6 % This determines the critical Mach number

%%-----------------------------------------------------------------------------------
%% Drag Estimites

%% Aircraft Parameters
W = 28229; % lbf
S_ref = 300; % ft^2
c = -0.1289; % Coefficient for fighter aircraft, given for S_wetrest equation,
% provided by Roskam's Aircraft Design Volume 1 (1985), Table 3.5.
d = 0.7506; % Coefficient for fighter aicraft, given for S_wetrest equation,
% provided by Roskam's Aircraf Design Volume 1 (1985), Table 3.5.
Cf = 0.0035; % Skin friction coefficient

M_lower = 0.1; % Lower Mach number
M_upper = 3.0; % Upper Mach number
% Generate Mach array
n = 100; % number of data points desired
mach_array = linspace(M_lower, M_upper, n);
mach_array = mach_array';
% Initializing some empty arrays
CD = zeros(length(mach_array), 1);
CD0_sub = zeros(length(CD), 1);
CD0_sup = zeros(length(CD), 1);
CD_sub = zeros(length(CD), 1);
CD_sup = zeros(length(CD), 1);
CD_wave = zeros(length(CD), 1);
CD_min = zeros(length(CD), 1);
k1_sup = zeros(length(CD), 1);

% Initializing some empty arrays
CD_mach = zeros(length(mach_array), 1);
CD0_sub = zeros(length(mach_array), 1);
CD0_sup = zeros(length(mach_array), 1);
CD_sub = zeros(length(mach_array), 1);
CD_sup = zeros(length(mach_array), 1);
CD_wave = zeros(length(mach_array), 1);
CD_min = zeros(length(mach_array), 1);
CL_mach = zeros(length(mach_array), 1);
k1_sup = zeros(length(mach_array), 1);
q = zeros(length(mach_array), 1);
V = zeros(length(mach_array), 1);


%%----------------------------------------------------------------------------------------------------------------------------
%% Compute S_wet
c = -0.1289; % Coefficient for fighter aircraft, given for S_wetrest equation, provided by Roskam's Aircraft Design Volume 1 (1985), Table 3.5.
d = 0.7506; % Coefficient for fighter aicraft, given for S_wetrest equation, provided by Roskam's Aircraf Design Volume 1 (1985), Table 3.5.
S_wet = 10^(c) * W^(d); % ft^2
CD0 = Cf * S_wet/S_ref;

% Generate velocity array
% Get local speed of sound at altitude
[T, rho, mu] = get_standard_atmosphere(h_ft); %GET ATMOSPHEREIC CONDITIONS @ ALTITUDE

% rho = slugs/ft^3
% mu = slugs/(ft*s)
% T = Rankine
% Get local speed of sound at 35000 ft
[~, a, ~, ~, ~, ~] = atmosisa(h_ft*0.3048); %GET SPEED OF SOUND

a = a/0.3048; % Convert m/s into ft/s
% Compute velocity for mach range
V = mach_array.*a; % ft/sec
% Compute dynamic pressure
q = 0.5 .* rho .* V.^2; % lbf/ft^2

% Compute e_osw
e_osw = (4.6*(1 - 0.033 * AR^(0.53)) * (cos(Lambda_LE*pi/180))^(0.1) - 3.3); 
% Equation from Brandt's work, "Aero" sheet, around cell I12.

% Get k1 subsonic
k1_sub = 1/(pi*e_osw*AR);

% Compute c_l_alpha
c_l_alpha = 0.1; % From Brandt's work, "Aero", ~cell E14)

% Compute e_notoswald
e_notoswald = (2)/(2 - AR + sqrt(4 + (AR^2) * (1 + tan(Lambda_LE*pi/180)^2) ) );

% Compute CL_alpha
CL_alpha = (c_l_alpha) / ( 1 + ( (57.3 * c_l_alpha) / (pi * e_notoswald * AR) ) )

% Get CL_minD
CL_minD = 0.027156 % From Brandt, "Aero" sheet, cell G20, hardcoded value
% CL_minD = CL_alpha * ( sqrt(NACA/1000) /2);
alpha_L0 = -1.049;
CL_minD = 0.054312 * -1 * alpha_L0/2

% Compute cl array
cl_airfoil = CL_alpha.*alpha;

% Compute k2 subsonic
k2_sub = -2 * k1_sub * CL_minD

% Compute CDmin
CD_min = Cf*S_wet/S_ref;

% Calculate CD0 for subsonic case
CD0_sub = CD_min + k1_sub.*CL_minD.^2;
CD0 = CD0_sub

% Compute k1 supersonic
k1_sup = ( (AR.*(mach_array.^2 - 1)./(4*AR*sqrt(mach_array.^2 - 1) - 2) ) .* cos(Lambda_LE*pi/180) );
k1_sup(imag(k1_sup)~=0)=0
% Compute k2 supersonic
k2_sup = 0
% Compute CD_wave
CD_wave = ( (4.5.*pi)./(S_ref) .* ((Amax - 0)/ Length)^2 .* E_wd .* (0.74 + 0.37 .* cos(Lambda_LE*pi/180)) .* (1 - 0.3 .* sqrt(mach_array - M_CD0_max)) );
CD_wave(imag(CD_wave)~=0.0)=0;
% Compute CD0 for supersonic case (which technically covers the entire range)
CD0_sup = CD0_sub + CD_wave;


%% Compute CL
% Assume lift = weight
CL_mach = W ./ (q * S_ref); % Brandt uses this form in their work. Dimensionally, it's lbf/ft^2 / lfb/ft^2, which all cancels out nicely
CL_angle = CL_alpha.*alpha;

%% Compute CD for level case
% Subsonic
CD_sub = CD0_sub + k1_sub.*CL_mach.^2 + k2_sub.*CL_mach;

% Supersonic
CD_sup = CD0_sup + k1_sup.*CL_mach.^2 + k2_sup.*CL_mach;
% Trim values that we don't need
CD_sub = CD_sub(mach_array<M_CD0_max); % Removes all values ABOVE MCD0_max because this isn SUBSONIC ONLY
CD_sub = CD_sub(CD_sub>0.0); % Removes all zero elements.

% Trim imaginary numbers
CD_sup(imag(CD_sup)~=0.0)=0; % Sets all imaginary components to 0.
CD_sup = CD_sup(mach_array>=M_CD0_max); % Remove all values below MCD0_max because
CD_sup = CD_sup(CD_sup>0.0); % Removes all zero elements. this is SUPERSONIC ONLY
CD_mach = [CD_sub; CD_sup] % Assemble CD array

% Compute CL/CD
CL_CD_mach = CL_mach./CD_mach;

%% Compute CD for various alpha ranges
% Subsonic
CD_angle = CD0_sub + k1_sub.*CL_angle.^2 + k2_sub.*CL_angle;

% Supersonic
% CD_sup = CD0_sub + k1_sup.*CL_alpha.^2 + k2_sup.*CL_alpha;
% Trim values that we don't need
% CD_sub = CD_sub(mach_array<M_CD0_max); % Removes all values ABOVE MCD0_max because this isn SUBSONIC ONLY
% CD_sub = CD_sub(CD_sub>0.0); % Removes all zero elements.
% Trim imaginary numbers
% CD_sup(imag(CD_sup)~=0.0)=0; % Sets all imaginary components to 0.
% CD_sup = CD_sup(CD_sup>0.0); % Removes all zero elements.
% CD_sup = CD_sup(mach_array>=M_CD0_max); % Remove all values below MCD0_max because this is SUPERSONIC ONLY
% CD_angle = [CD_sub]; % Assemble CD array

% Compute CL/CD
CL_CD_angle = CL_angle./CD_angle;

% Form complete k1, k2
index = (mach_array<M_Crit);
k1(index) = k1_sub;
k2(index) = k2_sub;
index = (mach_array>=M_Crit);
k1(index) = k1_sup(index);
k2(index) = k2_sup;

% Form complete CD0
CD0 = CD0_sup % This is because CD0_sup is negliigibly different from CD0_sub because CD0_sub is incorporated into the results

%% Generate plots
figure
hold on
grid on
plot(CL_mach, CD_mach, 'LineWidth', 2)
xlabel("C_L")
ylabel("C_D")
hold off
figure
hold on
grid on
title("C_L/C_D vs Mach Number (level) ")
plot(mach_array, CL_CD_mach, 'LineWidth', 2)
xlabel("Mach Number")
ylabel("C_L/C_D")
hold off
figure
hold on
grid on
title("C_L vs C_D (angle)")
plot(CL_angle, CD_angle, 'LineWidth', 2)
xlabel("C_L")
ylabel("C_D")
hold off

%% Generate table
% Index the Mach values close to the key ones
m_key = [0.1, 0.2, 0.5, M_CD0_max, M_Crit, 1.5, 2.0]; % Initialize array of key Mach values.

% Begin borrowed code:
% Unoriginal code. Source (top answer): https://www.mathworks.com/matlabcentral/answers/4249-round-towards-specific-values-in-an-array
for idx1=1:length(m_key);
    for idx2=1:length(mach_array);
        C(idx2,idx1)=m_key(idx1)-mach_array(idx2);
    end;
end
% Now find the index of the min values
[v,i]=min(abs(C));
% 'i' now contants the list of locations in B that corespond to the nearest
% A value
% mach_array(i)
% End borrowed code

% Form the table
T = table(mach_array(i), CD_mach(i), CD0(i), k1(i)', k2(i)');
T.Properties.VariableNames = {'Mach Number', 'CD0', 'CD', 'k1', 'k2'};
disp(T)

%% Compute CL for the range of values
% CL = W./(q.*S_ref);
CL = CL_alpha.*alpha;

% Uncambered wing
CD_uncambered = CD0 + k1_sub.*CL.^2

% Cambered wing
CD_cambered = CD0 + k1_sub.*(CL - CL_minD).^2

% Plot
figure
hold on
plot(CL, CD_uncambered, 'LineWidth', 2)
plot(CL, CD_cambered, 'LineWidth', 2)
title("CD vs CL")
xlabel("C_L")
ylabel("C_D")
legend("Uncambered", "Cambered")
grid on


%%-----------------------------------------------------------------------------------------------
%% Lift Computation

%% This script creates a .mat file that includes all of the problem givens

% Define mission segment data
mission = struct( 'Takeoff', struct('Altitude_ft', 0, 'Mach', 0.00, 'Time_min', [], 'Distance', [], 'DropPayload_lb', [], 'FixedPayload_lb', 700), ...
'Climb', struct('Altitude_ft', 40000, 'Mach', 0.87, 'Time_min', [], 'Distance', [], 'DropPayload_lb', [], 'FixedPayload_lb', 700), ...
'Cruise', struct('Altitude_ft', 40000, 'Mach', 0.87, 'Time_min', [], 'Distance', 190*6076.12, 'DropPayload_lb', [], 'FixedPayload_lb', 700), ...
'Dash', struct('Altitude_ft', 40000, 'Mach', 1.50, 'Time_min', [], 'Distance', 50*6076.12, 'DropPayload_lb', [], 'FixedPayload_lb', 700), ...
'Combat', struct('Altitude_ft', 25000, 'Mach', 0.87, 'Time_min', 2, 'Distance', [], 'DropPayload_lb', 4400, 'FixedPayload_lb', 700), ...
'Cruise2', struct('Altitude_ft', 40000, 'Mach', 0.87, 'Time_min', [], 'Distance', 250*6076.12, 'DropPayload_lb', [], 'FixedPayload_lb', 700), ...
'Loiter', struct('Altitude_ft', 10000, 'Mach', 0.30, 'Time_min', 20, 'Distance', [], 'DropPayload_lb', [], 'FixedPayload_lb', 700), ...
'Landing', struct('Altitude_ft', 0, 'Mach', [], 'Time_min', [], 'Distance', [], 'DropPayload_lb', [], 'FixedPayload_lb', 700));

% Define the propulsion struct
propulsion = struct( ...
'Takeoff', struct('TSFC', []), ...
'Climb', struct('TSFC', []), ...
'Cruise', struct('TSFC', 0.00019), ...
'Dash', struct('TSFC', 0.00061), ...
'Combat', struct('TSFC', 0.0004744), ...
'Cruise2', struct('TSFC', 0.00019), ...
'Loiter', struct('TSFC', 0.00019), ...
'Landing', struct('TSFC', []));
% 'Combat', struct('TSFC', 0.02847), ...;

% Define the aerodynamics struct
aero_mission = struct( ...
'Takeoff', struct('e', [], 'CD0', []), ...
'Climb', struct('e', [], 'CD0', []), ...
'Cruise', struct('e', 0.914, 'CD0', 0.027), ...
'Dash', struct('e', 0.914, 'CD0', 0.027), ...
'Combat', struct('e', 0.914, 'CD0', 0.017), ...
'Cruise2', struct('e', 0.914, 'CD0', 0.017), ...
'Loiter', struct('e', 0.914, 'CD0', 0.017), ...
'Landing', struct('e', [], 'CD0', []));

% Define an intermediate struct for air stuff
air = struct( ...
'Takeoff', struct('q', 0.00, 'a',1116.29), ...
'Climb', struct('q', 208.02,'a',967.94), ...
'Cruise', struct('q', 208.02,'a',967.94), ...
'Dash', struct('q', 618.37,'a',967.94), ...
'Combat', struct('q', 415.60,'a',1016.02), ...
'Cruise2', struct('q', 208.02,'a',967.94), ...
'Loiter', struct('q', 91.63,'a',1077.30), ...
'Landing', struct('q', 0.00,'a',1116.29));

% Define constraints
constraints = struct( ...
'MxMach', struct('W_Wto', 0.899666963, 'altitude', 36000, 'Mach', 1.6, 'n', 1, 'AB', 100, 'Ps', 0, 'CDx', 0), ...
'Cruise', struct('W_Wto', 0.899666963, 'altitude', 36000, 'Mach', 0.87,'n', 1, 'AB', 0, 'Ps', 0, 'CDx', 0), ...
'MaxAlt', struct('W_Wto', 0.899666963, 'altitude', 50000, 'Mach', 0.87,'n', 1, 'AB', 100, 'Ps', 0, 'CDx', 0), ...
'CmbtTm1',struct('W_Wto', 0.899666963, 'altitude', 20000, 'Mach', 0.87,'n', 4.5, 'AB', 100, 'Ps', 0, 'CDx', 0), ...
'CmbtTm2',struct('W_Wto', 0.899666963, 'altitude', 36000, 'Mach', 1.4, 'n', 1.4, 'AB', 100, 'Ps', 0, 'CDx', 0), ...
'Ps', struct('W_Wto', 0.899666963, 'altitude', 10000, 'Mach', 0.87,'n', 1, 'AB', 100, 'Ps', 500, 'CDx', 0));

% Define aero data
aero_constraints = struct( ...
'MxMach', struct('T_degR', 390.53, 'rho', 0.000707, 'CDo', 0.039317, 'k1', 0.213727, 'k2', 0, 'a', 968.61, 'V', 1549.78, 'q', 849.43), ...
'Cruise', struct('T_degR', 390.53, 'rho', 0.000707, 'CDo', 0.016996, 'k1', 0.116031, 'k2', -0.0063, 'a', 968.61, 'V', 842.69, 'q', 251.15), ...
'MaxAlt', struct('T_degR', 389.99, 'rho', 0.000363, 'CDo', 0.016996, 'k1', 0.116031, 'k2', -0.0063, 'a', 967.94, 'V', 842.11, 'q', 128.57), ...
'CmbtTm1',struct('T_degR', 447.49, 'rho', 0.001265, 'CDo', 0.016996, 'k1', 0.116031, 'k2', -0.0063, 'a', 1036.85,'V', 902.06, 'q', 514.75), ...
'CmbtTm2',struct('T_degR', 390.53, 'rho', 0.000707, 'CDo', 0.040614, 'k1', 0.219406, 'k2', -0.001, 'a', 968.61, 'V', 1356.06, 'q', 650.35), ...
'Ps', struct('T_degR', 483.09, 'rho', 0.001754, 'CDo', 0.016996, 'k1', 0.116031, 'k2', -0.0063, 'a', 1077.30,'V', 937.25, 'q', 770.59));

% Define thrust data (lapse ratios)
thrust = struct( ...
'MxMach', struct('dry_lapse', 0.298293, 'AB_lapse', 0.57698, 'throttle_lapse', 0.5770), ...
'Cruise', struct('dry_lapse', 0.27111, 'AB_lapse', 0.332642, 'throttle_lapse', 0.1711), ...
'MaxAlt', struct('dry_lapse', 0.138789, 'AB_lapse', 0.17029, 'throttle_lapse', 0.1703), ...
'CmbtTm1',struct('dry_lapse', 0.555662, 'AB_lapse', 0.681777, 'throttle_lapse', 0.6818), ...
'CmbtTm2',struct('dry_lapse', 0.357862, 'AB_lapse', 0.56558, 'throttle_lapse', 0.5566), ...
'Ps', struct('dry_lapse', 0.702727, 'AB_lapse', 0.85355, 'throttle_lapse', 0.8536));

% Define takeoff and landing data
TO = struct( ...
'Takeoff', struct('W_Wto', 1, 'V_Vstall', 1.2, 'mu', 0.03, 'Distance', 4000, 'air_density', 0.002377, 'CD0', 0.051996, 'AB_thrustlapse', 0.955053, 'CLmax', 1.27567));
Landing = struct( ...
'Landing', struct('W_Wto', 1, 'V_Vstall', 1.3, 'mu', 0.5, 'Distance', 4000, 'air_density', 0.002377, 'CD0', 0.061996, 'CLmax', 1.42591));

% Constants
AR = 3;
W_fixed = 5100;
save('mission_inputs.mat', 'mission', 'propulsion', 'aero_mission', 'aero_constraints', 'air', 'constraints', 'thrust', 'TO', 'Landing', 'W_fixed', 'AR');

% Initializing some variables
% Wing loading (from HW1 reference solution excel sheet)
mission.Takeoff.W_S = 104.59; % lbf/ft^2
mission.Climb.W_S = 99.36; % lbf/ft^2
mission.Cruise1.W_S = 97.20; % lbf/ft^2
mission.Dash.W_S = 94.40; % lbf/ft^2
mission.Combat.W_S = 92.09; % lbf/ft^2
mission.Cruise2.W_S = 80.55; % lbf/ft^2
mission.Loiter.W_S = 78.13; % lbf/ft^2
mission.Landing.W_S = 76.05; % lbf/ft^2
% Dynamic pressure (from HW1 reference solution Excel sheet)
mission.Takeoff.q = 210.70;
mission.Climb.q = 208.02; % lbf/ft^2
mission.Cruise1.q = 208.02; % lbf/ft^2
mission.Dash.q = 618.37; % lbf/ft^2
mission.Combat.q = 415.60; % lbf/ft^2
mission.Cruise2.q = 208.02; % lbf/ft^2
mission.Loiter.q = 91.63; % lbf/ft^2
mission.Landing.q = 76.16; % lbf/ft^2

%% ----------------------------
%% Design lift coefficient
% Obtain the lift coefficient for the mission segment for which you're
% designing your aircraft around. This segment should be the crux of
% everything, the most important segment to fulfil.
% Consider multiple mission segments...
% Takeoff
cl_takeoff = mission.Takeoff.W_S/mission.Takeoff.q
% Climb
cl_climb = mission.Climb.W_S/mission.Climb.q
% Cruise 1
cl_cruise1 = mission.Cruise1.W_S/mission.Cruise1.q
% Dash
cl_dash = mission.Dash.W_S/mission.Dash.q
% Combat
cl_combat = mission.Combat.W_S/mission.Combat.q
% Cruise 2
cl_cruise2 = mission.Cruise2.W_S/mission.Cruise2.q
% Loiter
cl_loiter = mission.Loiter.W_S/mission.Loiter.q
% Landing
cl_landing = mission.Landing.W_S/mission.Landing.q

%%---------------------------------------------------------------------------
%% Initializing some stuff
M_lower = 0.0; % Lower Mach limit
M_upper = 3.0; % Upper Mach limit
n = 100; % Number of data points
M = linspace(M_lower, M_upper, n); % Generate Mach array
M = M'; % Rotate so it's easier to view in the workspace variable list
C_l_alpha = 0.1; % From Brandt, "Aero" sheet, cell E14 I think.
d = 3.50 + 3.50; % Fuselage diameter, ft (taken from Brandt, "Geom," cells M89 & N89. Using maximum fuselage diameter.
b = 30; % Wingspan (ft)
A = 3.0; % Aspect ratio
lambda_LE = 40*pi/180; % Leading edge sweep (deg -> rad conversion)
lambda_TE = 0;
lambda_max_t = (lambda_LE - lambda_TE)/2;
S_ref_exposed = 196.23; % Exposed wing area (ft^2) (S_ref - wing_area_covered_by_fuselage)
S_ref = 300; % Wing reference area (ft^2)

%% Subsonic lift-curve slope
% Compute CL_alpha_subsonic for the range of Mach values
CL_alpha_subsonic = (2.*pi)./(sqrt(1 - M.^2));
% Clean up any imaginary values
CL_alpha_subsonic(imag(CL_alpha_subsonic)~=0)=0;
% Set any values in the transonic region to zero
CL_alpha_subsonic((0.85<M) & (M<1.2))=0;
% Plot
figure
plot(M, CL_alpha_subsonic, 'LineWidth', 2)
grid on
xlabel("Mach number")
ylabel("C_{L_{\alpha}}")
title("C_{L_{\alpha}} subsonic vs Mach Number")
%% Supersonic lift-curve slope
M_boundary = 1/cos(lambda_LE); % This computes the lowest Mach number considered for the supersonic case.
CL_alpha_supersonic = (4./sqrt(M.^2 - 1));

% Clean up any imaginary values
CL_alpha_supersonic(imag(CL_alpha_supersonic)~=0)=0;

% Set any values in the transonic region to zero
CL_alpha_supersonic((0.85<M) & (M<1.2))=0;

% Plot
figure
plot(M, CL_alpha_supersonic, 'LineWidth', 2)
grid on
xlabel("Mach number")
ylabel("C_{L_{\alpha}}")
title("C_{L_{\alpha}} supersonic vs Mach Number")

%% Using general form thing I guess?
Beta = sqrt(1 - M.^2);
eta = C_l_alpha./((2.*pi)./Beta);
F = 1.07*(1 + d/b)^2;
CL_alpha = (2.*pi.*A)./(2 + sqrt(4 + ((A.^2 .* Beta.^2)./(eta.^2)) .* (1 + (tan(lambda_max_t)^2)./(Beta.^2)))).*(S_ref_exposed./S_ref).*(F);
% Accurate up to the drag divergence Mach number for the swept wing

% Clean up any imaginary values
CL_alpha(imag(CL_alpha)~=0)=0;

% Set any values in the transonic region to zero
CL_alpha((0.85<M) & (M<1.2))=0;

% Plot
figure
plot(M, CL_alpha, 'LineWidth', 2)
grid on
xlabel("Mach number")
ylabel("C_{L_{\alpha}}")
title("C_{L_{\alpha}} general vs Mach Number")

%% Interpolation
% Interpolate between Mach 0.85 & Mach 1.2.
x = find(CL_alpha ~= 0);
y = CL_alpha(x);
CL_alpha = interp1(x, y, 1:numel(CL_alpha), 'spline', 'extrap');

% Plot
figure
plot(M, CL_alpha, 'LineWidth', 2)
grid on
xlabel("Mach number")
ylabel("C_{L_{\alpha}}")
title("C_{L_{\alpha}} interpolated vs Mach Number")

%%--------------------------------------------------------
% Initializing some constants
cl_alpha = 0.1;
lambda_tmax = 40/2; % Degrees
AR = 3.0;
AR_t = 3.0;
S_ref = 300; % ft^2
S_pitch = 108; % ft^2
S_strakes_exposed = 20; % ft^2
wing_root_chord = 16.29; % ft
wing_tip_chord = 3.71; % ft
pitch_xloc = 36.0; % ft
pitch_tipchord= 2.22; % ft
pitch_rootchord = 9.78; % ft
wings_xloc = 17.79; % ft
wings_tipchord = 3.71; % ft
wings_rootchord = 16.29; % ft
lambda = 0.23;
z_h = 0.0;
b = 30.0;
c_avg = (wing_root_chord + wing_tip_chord)/2

wing_c = wings_xloc + 0.25*wings_rootchord
pitch_c = pitch_xloc + 0.25*pitch_rootchord;

l_h = (pitch_c - wing_c)

% Compute e
e_notoswald = 2/(2 - AR + sqrt(4 + AR^2 * ( 1 + tand(lambda_tmax)^2)))

% Compute CL_alpha main wing
CL_alpha = cl_alpha/(1 + (57.3*cl_alpha)/(pi * e_notoswald * AR))

% Compute CL_alpha_t pitch trim surface
CL_alpha_t = cl_alpha/(1 + (57.3*cl_alpha)/(pi * e_notoswald * AR_t))

% Compute deltaepsilon_deltaalpha
deltaepsilon_deltaalpha = (21 * CL_alpha)/(AR^(0.5)) * (c_avg/l_h)^(0.25) * ((10 - 3*lambda)/7) * (1 - z_h/b)

% Compute CL_alpha for whole aircraft WITH strakes
CL_alpha_whole_strakes = CL_alpha * (S_ref + S_strakes_exposed) / S_ref

% Compute CL_alpha for whole aircraft
CL_alpha_whole = CL_alpha_whole_strakes + CL_alpha_t*(1 - deltaepsilon_deltaalpha) * (S_pitch/S_ref)

%-------------------------------------------------------------------------------
%% Computing various CLmax values for the F-16A Falcon.
% Initializing some constants

NACA = 1404; % Given NACA 4-digit airfoil number that doesn't really exist buuut gives uuus the right valuues anyway
alpha_max = 15; % Maximum AOA assumed (deg)
% alpha_L0 = NACA/1000; % Zero-lift AOA (deg)
alpha_L0 = -1.0469; % Zero-lift AOA (deg) (obtained from Xfoil, inviscid, M 0.2)
cl_alpha = 0.1;
lambda_tmax = 40/2; % Degrees
AR = 3.0;
AR_t = 3.0;
S_ref = 300; % ft^2
S_pitch = 108; % ft^2
S_strakes = 20; % ft^2
wing_root_chord = 16.29; % ft
wing_tip_chord = 3.71; % ft
c_avg = (wing_root_chord + wing_tip_chord)/2;
pitch_xloc = 36.0; % ft
pitch_tipchord= 2.22; % ft
pitch_rootchord = 9.78; % ft
wings_xloc = 17.79; % ft
wings_tipchord = 3.71; % ft
wings_rootchord = 16.29; % ft
lambda = 0.23;
z_h = 0.0;
b = 30.0;

wing_tip_quarter_chord = 15.0; % ft
TE_flap_top_view_wingtip_LE_left_y = 11.25; % ft
wing_tip_chord = 3.71; % ft
wing_root_chord = 16.29; % ft
wing_wingtip_LE = 15.0; % ft
TE_flap_left_y = 3.50; % ft
wing_exposed_root_chord = 13.36; % ft
CL_alpha = 0.054312;
CL_alpha_whole = 0.0615; % CL_alpha for the entire aircraft + strakes
lambda_TE = 0;
l_h = (pitch_xloc - 0.25*wings_tipchord - wings_xloc - 0.25*wings_rootchord);

% Compute alpha_a_max
alpha_a_max = alpha_max - alpha_L0;

% Compute CLmax clean
CL_max_clean = CL_alpha_whole * alpha_a_max
% Compute area of horizontal stabilizer
S_h = 108; % ft^2

% Compute Delta CLmax (max trimmable Delta CL)
Delta_CL_max = S_h/S_ref * l_h / (0.5 * wing_root_chord);

% Compute the flapped area of the main wings
S_flapped = (S_ref - ( wing_tip_quarter_chord - TE_flap_top_view_wingtip_LE_left_y ) * 2 * ...
    ( wing_tip_chord + ( wing_root_chord - wing_tip_chord)*(wing_tip_quarter_chord - ...
    TE_flap_top_view_wingtip_LE_left_y)/(wing_wingtip_LE)) - (TE_flap_left_y*(wing_root_chord + wing_exposed_root_chord)))

% Compute Delta alpha for landing flaps
Delta_alpha_landing = 15*S_flapped / (S_ref * cos(lambda_TE*pi/180));

% Compute the Delta CLmax for landing
Delta_CL_max_landing = CL_alpha_whole * Delta_alpha_landing;

% Compute CLmax during landing
CL_max_landing = CL_max_clean + Delta_CL_max_landing

% Compute CLmax during takeoff
CL_max_takeoff = CL_max_clean + (2/3)*Delta_CL_max_landing

%%----------------------------------------------------------------------------------------------------
%% Estimate the CL_max for various situations. Clean, according to Raymer,
%% Initializing some constants
lambda_LE = 40; % Leading edge sweep (deg)
lambda = 0.23; % Taper ratio
AR = 3.0; % Aspect ratio
NACA = 1404; % Equivalent NACA 4-digit

% Compute the main wing's sweep at the quarter chord
Lambda_quarterchord = atan(tan(lambda_LE*pi/180) - ( 1 - lambda)/(AR*(1+lambda)))*180/pi % Compute the quarter-chord sweep angle (deg)

% Compute the 2-D Cl_max for the loiter, cruise, and dash segments
% Loiter
mission.Loiter.W_S = 78.13; % lbf/ft^2
mission.Loiter.q = 91.63; % lbf/ft^2
cl_max_loiter = mission.Loiter.W_S/mission.Loiter.q

% Cruise 1
mission.Cruise1.W_S = 97.20; % lbf/ft^2
mission.Cruise1.q = 208.02; % lbf/ft^2
cl_max_cruise = mission.Cruise1.W_S/mission.Cruise1.q

% Dash
mission.Dash.W_S = 94.40; % lbf/ft^2
mission.Dash.q = 618.37; % lbf/ft^2
cl_max_dash = mission.Dash.W_S/mission.Dash.q

%% Compute the CL_max for the aircraft for various segments
% Loiter
CL_max_clean_loiter = 0.9*cl_max_loiter * cosd(Lambda_quarterchord); % Eq 12.15 (valid for subsonic aircraft of "moderate sweep").

% Cruise
CL_max_clean_cruise = 0.9*cl_max_cruise * cosd(Lambda_quarterchord)

% Dash
CL_max_clean_dash = 0.9*cl_max_dash * cosd(Lambda_quarterchord);

% Note: "Moderate sweep" isn't defined. Therefore, I'll assume it's < 30
% degrees. Therefore, the F-16 is not moderately swept; these equations
% won't accurately estimate the CL_max_clean.
%% Compute CL_max for high aspect ratio ( M = 0.2 )
% Set AR = 5.0
AR = 5.0;
% Determine the sharpness parameter
% Should be tabulation.
% NACA 4 digit
% t_c = (NACA/1000)/;
t_c = 0.04
Delta_y = 21.3 * t_c % You should use the NACA 6-digit

% Getting CL_max/cl_max
% Tabulating...
% Lambda_LE = 40 degrees
CL_max_cl_max = 1.1; % Tabulation; give generous error margin

% Compute CL_max for the airfoil for Machs 0.2, 0.4, & 0.6
% Obtain Cl_max for Mach 0.2
Cl_max_m02 = 0.1;

% Obtain Delta_CL_max for Mach 0.2, 0.4, 0.6, from Fig 12.9
Delta_CL_max = [-0.0001; -0.025; -0.05];

% Compute CL_max_highAR for Mach numbers
CL_max_highAR = Cl_max_m02.*(CL_max_cl_max) + Delta_CL_max % According to Roskam, I should be getting values between 1.2 - 1.8.

% This is to demonstrate the results of using an incorrect approach to
% determining CL_max

%% Compute alpha_CL_max for high AR wing
% Obtain CL_max_cl_max ratio (done)

% Obtain alpha_L0
% alpha_L0 = -NACA/1000;
alpha_L0 = -1.0469; % From Xfoil, inviscid solution, Mach 0.2.

% Obtain Delta_alpha_CL_max from figure 12.10
% Tabulating...
Delta_alpha_CL_max = 7.25; % Tabulation; give generous error margin

% Obtain CL_alpha for M 0.2, 0.4, 0.6
% CL_alpha = [0.0801; 0.0780; 0.0727];
CL_alpha = 0.054312; % From Brandt's work

% Computing alpha_CL_max for high AR wing for Machs 0.2, 0.4, & 0.6
alpha_CL_max_highAR = CL_max_highAR./CL_alpha + alpha_L0 + Delta_alpha_CL_max

%% Low AR wing
% Determine if wing is low AR
AR = 3.0; % Set AR to F-16A's true value.

% Obtain C1 & C2
C1 = 0.5;
C2 = 0.6;

% Compute AR_test
AR_Limit_Test = 3/( (C1 + 1) * (cos(lambda_LE*pi/180)));
disp("AR limit: " + AR_Limit_Test)
if AR <= AR_Limit_Test
    disp("Current AR (" + AR + ") below AR limit. Use low-AR equations!")
elseif AR > AR_Limit_Test
    disp("Current AR (" + AR + ") above AR limit. Use high-AR equations!")
end

%%---------------------------------------------------------------------------------
%% Estimating Delta_CL_max under the influence of various types of flaps.
% Initializing some constants
% Main wing
wing_c_root_exposed = 13.36; % Exposed root chord length (ft)
wing_c_root = 16.29; % Geometric root chord length (ft)
wing_c_mean = (wing_c_root_exposed + wing_c_root)/2; % Using the mean to account for students using the tip and/or root chord lengths

% TE flaps
flap_c_root = 2.02; % ft
flap_c_tip = 1.08; % ft
flap_c_mean = (flap_c_root + flap_c_tip)/2; % Using the mean to account for students using the tip and/or root chord lengths
alpha_0L = -1.0469; % Zero-lift angle of attack for the 2-D airfoil (deg)
S_ref = 300; % Reference planform area (ft^2)
CL_max_clean = 0.9869; % Assumed value (Brandt)
Delta_alpha_L0_airfoil = 0; % Symmetric airfoil
%% Computing Sflapped
% Need width of aileron/elevon
% Rectangle A...
A_w = 11.25 - 3.50;
A_h = 34.079 - 30.37;
A_Area = A_w*A_h;
% Triangle B...
B_w = A_w;
B_h = wing_c_root - A_h;
B_Area = (0.5)*B_w*B_h;
% Computing area
S_flapped = A_Area + B_Area
% Note: For some reason, these are the exact same dimensions as the
% aileron/elevon in Brandt's workbook ("Geom" page, TE Flap top view and
% aileron/elevon work).
%% Computing hinge line's sweep angle
% Forming triangle C from triangular portion of aileron/elevon...
A_h = 34.079 - 32.079
C_h = A_h - (34.08 - 33.0)
Lambda_HL = atan(C_h/A_w)
%% Flaps
% Plain
Delta_Cl_max_plain = 0.9;
% Split
Delta_Cl_max_split = 0.9;
% Slotted flap
Delta_Cl_max_slotted = 1.3;
% Slotted fowler flap
c_extended = wing_c_root + 2.02;
Delta_Cl_max_fowler = 1.3 * c_extended/wing_c_root;
% Double slotted flap
Delta_Cl_max_doubleslottedflap = 1.6 * c_extended/wing_c_root;
% Triple slotted flap
Delta_Cl_max_tripleslottedflap = 1.9 * c_extended/wing_c_root;
%% Computing Delta_CL_max for each flap
% Plain
Delta_CL_max_plain = Delta_Cl_max_plain*(S_flapped/S_ref)*cos(Lambda_HL);
% Split
Delta_CL_max_split = Delta_Cl_max_split*(S_flapped/S_ref)*cos(Lambda_HL);
% Slotted
Delta_CL_max_slotted = Delta_Cl_max_slotted*(S_flapped/S_ref)*cos(Lambda_HL);
% Fowler
Delta_CL_max_fowler = Delta_Cl_max_fowler*(S_flapped/S_ref)*cos(Lambda_HL);
% Double slotted
Delta_CL_max_doubleslotted = Delta_Cl_max_doubleslottedflap*(S_flapped/S_ref)*cos(Lambda_HL);
% Triople slotted
Delta_CL_max_tripleslotted = Delta_Cl_max_tripleslottedflap*(S_flapped/S_ref)*cos(Lambda_HL);
%% Compute the new CL maxes for each flap.
% Plain
CL_max_plain = CL_max_clean + Delta_CL_max_plain;
% Split
CL_max_split = CL_max_clean + Delta_CL_max_split;
% Slotted
CL_max_slotted = CL_max_clean + Delta_CL_max_slotted;
% Fowler
CL_max_fowler = CL_max_clean + Delta_CL_max_fowler;
% Double slotted
CL_max_doubleslotted = CL_max_clean + Delta_CL_max_doubleslotted;
% Triple slotted
CL_max_tripleslotted = CL_max_clean + Delta_CL_max_tripleslotted;
% For the students, consider the following:
% Which flap gives the highest CL_max?
% What mission segment would this be best suited for?
% Consider the flaps that the F-16 actually uses. Are they the same as what
% gives the highest values, here? Why do you think they went with something
% else?

%%--------------------------------------------------------------------------------------
% Rstimating Delta_CL_max under the influence of various types of flaps.
% Initializing some constants
% Main wing
wing_c_root_exposed = 13.36; % Exposed root chord length (ft)
wing_c_root = 16.29; % Geometric root chord length (ft)
wing_c_mean = (wing_c_root_exposed + wing_c_root)/2; % Using the mean to account for students using the tip and/or root chord lengths

% TE flaps
flap_c_root = 2.02; % ft
flap_c_tip = 1.08; % ft
flap_c_mean = (flap_c_root + flap_c_tip)/2; % Using the mean to account for students using the tip and/or root chord lengths
alpha_0L = -1.0469; % Zero-lift angle of attack for the 2-D airfoil (deg)
S_ref = 300; % Reference planform area (ft^2)
CL_max_clean = 0.9869; % Assumed value (Brandt)

%% Computing Sflapped
% Computing area
S_flapped = S_ref; % This is because the LE flaps extend the ENTIRE LENGTH OF THE WING

%% Computing hinge line's sweep angle
% Assuming quarter chord sweep angle for simplicity's sake
Lambda_HL = 32.228 * pi/180; % deg to rad

%% Computing Delta_Cl_max for various TE devices
% Fixed slot
Delta_Cl_max_fixedslot = 0.2;

% Leading edge flap
Delta_Cl_max_leadingedgeflap = 0.3;

% Kruger flap
Delta_Cl_max_krugerflap = 0.3;

% Slat
c_extended = wing_c_root + 2.02;
Delta_Cl_max_slat = 0.4 * c_extended/wing_c_root;

%% Computing Delta_CL_max for various LE devices
% Fixed slot
Delta_CL_max_fixedslot = Delta_Cl_max_fixedslot*(S_flapped/S_ref) * cos(Lambda_HL);

% Leading edge flap
Delta_CL_max_leadingedgeflap = Delta_Cl_max_leadingedgeflap*(S_flapped/S_ref) * cos(Lambda_HL);

% Kruger flap
Delta_CL_max_krugerflap = Delta_Cl_max_krugerflap*(S_flapped/S_ref) * cos(Lambda_HL);

% Slat
Delta_CL_max_slat = Delta_Cl_max_slat*(S_flapped/S_ref) * cos(Lambda_HL);

%% Computing CL_max for various LE devices
% Fixed slot
CL_max_fixedslot = CL_max_clean + Delta_CL_max_fixedslot;

% Leading edge flap
CL_max_leadingedgeflap = CL_max_clean + Delta_CL_max_leadingedgeflap;

% Kruger flap
CL_max_krugerflap = CL_max_clean + Delta_CL_max_krugerflap;

% Slat
CL_max_slat = CL_max_clean + Delta_CL_max_slat;

%%--------------------------------------------------------------------------------
%% Sweep angle Mach Cone Angles
% Constants
lambda = 0.23;
M = 1.50;
AR = 3.0;

% Tabulating values for LE sweep angle
Lambda_LE_tabulated = 45 % Deg (give generous tolerance)

% Computing quarter-chord sweep angle
Lambda_quarterchord_tabulated = atan(tan(Lambda_LE_tabulated*pi/180) - ( (1 - lambda)/(AR * (1+lambda))))*180/pi

% Calculating value for LE sweep angle
Lambda_LE_calculated = 90 - asin(1/M)*180/pi

% Calculating value for quarter-chord sweep angle
Lambda_quarterchord_calculated = atan(tan(Lambda_LE_calculated*pi/180) - ( (1 - lambda)/(AR * (1 + lambda))))*180/pi

% Estimating Mach cone angle
Mach_cone_angle = asin(1/M)*180/pi