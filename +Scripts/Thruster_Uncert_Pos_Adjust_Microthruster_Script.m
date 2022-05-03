%% Sub-Microthruster Fusion Technique Simulation: Thruster Position Tolerancing
% Tests ranging uncertainty in the geometry of the spacecraft; leaves all
% other values and disturbances off except for measurement and process
% noise.
% Jonathan Hood, 5/3/2022
% Simulate an on-orbit sub-microthruster evaluation via the proposed Fusion
% Technique, including orbital dynamics and present perturbations.
%% Spacecraft Geometry
% To be modified and iterated to meet performance requirements (10^-7 or >
% thrust resolution)

clear
clc
close all

% Assume one thou for each directional tolerance
mech_tol = 1;
thrust_loc_arr = zeros([3 100]);
for kk = 1:100

% Set spacecraft tolerancing values
sc_tol = 1/1000*0;% spacecraft +/- tolerance in m

% w_b_ECI0 = [0.0; 0; 0]; % initial angular velocity rad/s ARBITRARY from detumble

w_b_ECI0 = [0.05; -0.05; 0.05]; % initial angular velocity rad/s
quat_b_ECI0 = [sin(1/2)/sqrt(3)*ones([3 1]); cos(1/2)]; % initial quaternion ARBITRARY from detumble
init_state = [w_b_ECI0; quat_b_ECI0];

units = 1;               % number of cubes in configuration
cube_max_mass = 1.33;    % Maximum mass of a single cube in kg
cube_len = 0.1;          % length of a side of a cube in m
m = cube_max_mass*units; % total allowable mass of s/c in kg
L = cube_len*units + randn(1)*sc_tol*mech_tol; % total length of the spacecraft 
                         % in meters
h = 0.1+ randn(1)*sc_tol*mech_tol;% height of a side of a cube in m
w = 0.1+ randn(1)*sc_tol*mech_tol;% width of a side of a cube in m

geom = [L h w]; % create geometry vector of Length (x), Height (z), Width (y)

% Length parallel to x-axis
% Width parallel to y-axis
% Height parallel to z-axis
% Inertia values are defined by the face the vector travels out of
% Ixx is then a function of y and z, or width and height
% Iyy is a function of x and z, length and height
% Izz is a function of x and y, length and width

J = 1/12*m*[w^2+h^2 0       0;
            0       h^2+L^2 0;
            0       0       w^2+L^2]; % s/c mass moment of inertia kg*m^2

% Longest dimension along x-axis
% Thruster placement is at edge of longest dimension, i.e at either +L/2 or
% -L/2

% Calculate constant torque vector
% thrust_loc_tol = 1/1000; % Set thruster location tolerance in in
thrust_tol = 1/1000*1;% thruster loc +/- tolerance in m

lever_arm_len = L/2; % lever arm length in m
theor_thrust = 1e-7; % predicted value of desired thrust to be measured
thruster_loc = [L/2+randn(1)*thrust_tol;thrust_tol*randn(1);thrust_tol*randn(1)]; % location of thruster at far edge of x-axis
thrust_loc_arr(:,kk) = thruster_loc;
az_tol = 0.001*0; % Error in thruster force vector in the xy plane
dec_tol = 0.001*0; % Error in thruster force vector about the z-plane
theta = 0; % optimal azimuthal theta is aligned with the y-axis
gamma = 0; % optimal declination gamma is 0, or in the xy plane
thrust_vec_unit = [sind(theta+az_tol); cosd(az_tol+theta); sind(dec_tol+gamma)];
thrust_vec_unit = thrust_vec_unit./norm(thrust_vec_unit);
thrust_vec = theor_thrust.*thrust_vec_unit; % Constant thrust vector in the positive y-dir
theor_torque = cross(thruster_loc,thrust_vec); % torque vector is cross of thruster position and thrust vector

frame_scale = 0.1;

figure
quiver3([0;0;0],[0;0;0],[0;0;0],[frame_scale; 0; 0],[0; frame_scale; 0],[0; 0; frame_scale],'LineWidth',4)
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
quiver3(0,0,0,thruster_loc(1),thruster_loc(2),thruster_loc(3))
quiver3(thruster_loc(1),thruster_loc(2),thruster_loc(3),thrust_vec(1).*1e7,thrust_vec(2).*1e7.*lever_arm_len,thrust_vec(3).*1e7,'-g')
quiver3(0,0,0,theor_torque(1).*1e7,theor_torque(2).*1e7,theor_torque(3).*1e7)
title('Thrust Vector in Body Frame')
legend('Body Frame','Thruster Position','Thrust Vector','Torque Vector')
view(136,23)
drawnow
exportgraphics(gcf, 'thruster_representation.png')
%% Orbital Motion Parameters
% Assumed to start in a rideshare SSO orbit.

% Calculate position
alt = 550;          % s/c altitude in km (from SpaceX average rideshare SSO altitude)
rE = 6378;          % radius of the Earth in km SPECIFIC
mu = 398600;        % std. grav param km^3/s^2
rmag = rE+alt;      % s/c position in km
inc = 97;           % relatively ARBITRARY SSO inclination in deg
ecc = 0;            % assume circular orbit
h = sqrt(rmag*mu);  % relative angular momentum of circular
RA = 0;             % let RAAN be 0
w = 0;              % let argument of perigee be 0
TA = 0;             % let TA be 0

coe = [h, ecc, RA, inc, w, TA];         % Create COES vector
[r,v] = Functions.sv_from_coe(coe,mu);  % get initial r & v vector in km and km/s
COES = Functions.COES(r,v,mu);          % Extract more detailed COES


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Initial Conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set desired date for simulation run and calculate Julian date
timeString1 = '2021-02-06 23:59:00'; % Time of start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create MATLAB datetime object in PST
t1 = datetime(timeString1,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss');
t1_PST = t1;

% Change t1 to UTC for math
t1.TimeZone = 'UTC';
JD_start = juliandate(t1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Initial Rotations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial b_ECI
C_b_ECI0 = Functions.quat2rot(quat_b_ECI0); % body ECI rotation matrix
euler_b_ECI0 = Functions.rot2euler(C_b_ECI0); % body ECI euler angles

% %%%%%%%%%%%%% CALCULATE LVLH ROTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
z = -r/norm(r);                                % z basis vector LVLH
y = -cross(r,v)/norm(cross(r,v));              % y basis vector LVLH
x = cross(y,z);                                % x basis vector LVLH

C_LVLH_ECI0 = [x;y;z];                             % get ECI to LVLH rotation matrix 
quat_LVLH_ECI0 = Functions.rot2quat(C_LVLH_ECI0);   % get ECI to LVLH quaternion
euler_LVLH_ECI0 = Functions.rot2euler(C_LVLH_ECI0); % get ECI to LVLH euler angles

% Initial b_LVLH
C_b_LVLH0 = C_b_ECI0*C_LVLH_ECI0';                     % get b to LVLH rotation matrix 
quat_b_LVLH0 = Functions.rot2quat(C_b_LVLH0);         % get b to LVLH quaternion
euler_b_LVLH0 = Functions.rot2euler(C_b_LVLH0);       % get b to LVLH euler angles

%% Set Torque Magnitude of Unmodeled Perturbations

dist_on = 0;

% From sim, know drag is essentially negligible (1e-26)
drag_unmodeled_unnormalized = randn([3 1]);
drag_unmodeled = 1e-26.*drag_unmodeled_unnormalized./norm(drag_unmodeled_unnormalized).*0;

% From sim, know SRP is essentially negligible (1e-24)
srp_unmodeled_unnormalized = randn([3 1]);
srp_unmodeled = 1e-24.*srp_unmodeled_unnormalized./norm(srp_unmodeled_unnormalized).*0;

% From sim, know Mag is  nonnegligible (1e-5); must model to at least 1e-8
mag_field_unmodeled_unnormalized = randn([3 1]);
mag_field_unmodeled = 1e-9.*mag_field_unmodeled_unnormalized./norm(mag_field_unmodeled_unnormalized).*dist_on;

% From sim, know GG is nonnegligible (1e-7); must model to at least 1e-8,
% adjust as needed
earth_and_oblation_mag = 2e-9.*dist_on;
% tidal_mag = 1e-9.*dist_on;
% sun_and_moon_mag = 1e-14.*dist_on;

GG_mag = earth_and_oblation_mag;
gg_unmodeled_unnormalized = randn([3 1]);
grav_grad_unmodeled = GG_mag.*gg_unmodeled_unnormalized./norm(gg_unmodeled_unnormalized).*dist_on;

% General relativity
lens_thirring_mag = 1e-18.*dist_on;
schwarzchild_mag = 1e-15.*dist_on;
de_sitter_mag = 1e-18.*dist_on;

GR_mag = lens_thirring_mag+schwarzchild_mag+de_sitter_mag;
gr_unmodeled_unnormalized = randn([3 1]);
gr_eff_unmodeled = GG_mag.*gr_unmodeled_unnormalized./norm(gr_unmodeled_unnormalized).*0;

% Onboard (unknown)
em_mag = 1e-9.*dist_on;
temp_mag = 1e-9.*dist_on;
connectors_mag = 1e-9.*dist_on;

on_mag = em_mag+temp_mag+connectors_mag;
onboard_unmodeled_unnormalized = randn([3 1]);
onboard_unmodeled = on_mag.*onboard_unmodeled_unnormalized./norm(onboard_unmodeled_unnormalized);

%% VSD UKF and ADD Parameters
% Section to set initial conditions, including covariance matrices and
% noise, for the variable-state dimension unscented Kalman filter.

% Modified from HW_5 for NCEKF test
P0 = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.1]); % initial covariance

% Measurement noise
sigma_m = 0.0000005; % std deviation magnetometer measurement noise
var_m = sigma_m^2; % magnetometer measurement noise variance
sigma_s = 0.005; % std deviation sun sensor measurement noise
var_s = sigma_s^2; % sun sensor measurement noise variance
sigma_g = 0.001; % gyroscope measurement noise
var_g = sigma_g^2; % gyroscope measurement noise variance
% R_k = diag([sigma_m^2,sigma_m^2,sigma_m^2,sigma_s^2,sigma_s^2,sigma_s^2]); % initial measurement noise covariance

% Test initial mag sensor covariance
RR_k = diag([sigma_m^2,sigma_m^2,sigma_m^2]);

% % Test initial sun sensor covariance
R_k = diag([sigma_s^2,sigma_s^2,sigma_s^2]);

% Initial gyroscope covariance
RRR_k = diag([var_g,var_g,var_g]);

% Process Noise
sigma_q = 0.001; % std deviation process noise
var_q = sigma_q^2; % process noise variance
rho = var_q;
Q_k = rho*eye(); % process noise covariance
% If rho >> norm(R), follow measurements closely
% If rho << norm(R), follow system model closely
% L_k = [delta_t*inv(J); zeros(3); zeros([1 3])];

% Gyroscope
g_sf_cc = [1 0 0; 0 1 0; 0 0 1];
g_bias = [0 0 0];
g_sen = [0 0 0];
wn = 190; % natural frequency
zn = 0.707; % damping ratio
g_sath = inf;
g_satl = -inf;
dtype_g = 0; % indicate include second order dynamics
%% Technique Specifics & Simulation Run
% Technique depends on precise measurement of spacecraft attitude, e.g.
% deflection of the opposing end of the lever arm.

% Assuming center of mass is in the geometric center of the spacecraft,
% lever arm is half the spacecraft length

dev_en = [1 0 1]; % Determines which attitude determination devices will be considered in the sim run
                  % Use 0 to disable, 1 to enable

% Turn thrust on/off
thrust_on = 1;
% Set measurement timesteps
sun_meas_frequ = 1;%frequ_range(jj); % sun sensor measurement frequency in Hz
mag_meas_frequ = 1; % magnetometer measurement frequ in Hz
gyro_meas_frequ = 50;%frequ_range(kk); % gyroscope measurement frequ in Hz
g_Ts = 1/gyro_meas_frequ; % gyroscope sample time for gyro functions
ts = g_Ts;

meas_frequ_arr = [sun_meas_frequ mag_meas_frequ gyro_meas_frequ];
% Set integration timestep; must be integer divisor of measurement
% frequencies
max_frequ = max(meas_frequ_arr);
delta_t = 1/max_frequ; % time step between measurements & prediction time step
% Set simulation time
simtime = 200;       % desired simulation runtime in sec
% Set initial state guess
w_hat0 = [0;0;0]; % initial angular velocity
quat_hat0 = [0;0;0;1]; % initial quaternion
x_hat0 = [w_hat0; quat_hat0]; % initial state

% Relative tolerance adjusted from default 1e-3 to 1e-12 to improve
% resolution
theor_torque = theor_torque.*thrust_on;
load_system("+Sims/attempt_2_Full_Microthruster_Simulation"); %load simulation
set_param('attempt_2_Full_Microthruster_Simulation', 'StopTime', 'simtime') %set sim to stop at desired runtime
outKF = sim('attempt_2_Full_Microthruster_Simulation');      % Run Simulink model

%% Plot MSD for all three systems
% load('outKF_Full_Sim_data.mat')
w_b_ECI = outKF.w_b_ECI.Data;
domega_b_ECI = outKF.domega_b_ECI.Data; % Extract differential angular velocity for torque calculations
x_hat_EKF = outKF.x_hat_EKF.Data;
x_hat_UKF = outKF.x_hat_UKF.Data;

% x_hat_k = outKF.x_hat_k.Data;
w_b_ECI_ekfilt = x_hat_EKF(:,1:3);
quat_b_ECI_ekfilt = x_hat_EKF(:,4:7);
w_b_ECI_ukfilt = x_hat_UKF(:,1:3);
quat_b_ECI_ukfilt = x_hat_UKF(:,4:7);
quat_b_ECI = outKF.quats_b_ECI.Data;
omega_thrust = outKF.omega_thrust.Data;
t = outKF.tout;
% t1= outKF.x_hat_EKF.Time;
t1= outKF.x_hat_UKF.Time;
sun_t = outKF.y_UKF_sun.Time;
gyro_t = outKF.y_UKF_gyro.Time;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRUST CALCS AND PLOTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% With Noise EKF Thrust Measurement

Tbb_EKF = zeros([3 length(t1)-1]);
Tbb_EKF_norms = zeros([1 length(t1)]);
for i = 2:length(t1)
    omega_b_b_2 = w_b_ECI_ekfilt(i,:)';
    omega_b_b_1 = w_b_ECI_ekfilt(i-1,:)';
    domega = (omega_b_b_2-omega_b_b_1)./delta_t; % Calculate change in omega, weighted by timestep
    Tbb_EKF(:,i) = J*domega+Functions.vec2cross(omega_b_b_1)*J*omega_b_b_1;
    
    quat_b_ECI_2 = quat_b_ECI_ekfilt(i,:)';
    C_b_ECI_2 = Functions.quat2rot(quat_b_ECI_2);
    Tbb_EKF(:,i) = C_b_ECI_2*Tbb_EKF(:,i);
    
    Tbb_EKF_norms(i) = norm(Tbb_EKF(:,i));
end

strIndex = 51; % at four seconds sim time to let filters converge

Tbb_EKF_calc = Tbb_EKF(:,strIndex:end); % Truncate torque values to approximately where convergence occurs
mean_torque_EKF = mean(Tbb_EKF_calc,2); % Find mean torque
mean_torque_uncert_EKF = std(Tbb_EKF_calc,0,2)./sqrt(length(Tbb_EKF_calc)); % Calculate uncertainty in each component
norm_uncert_EKF = norm(mean_torque_uncert_EKF); % Error in torque vector norm is never greater than the norm of the errors
thrust_uncert_EKF = norm_uncert_EKF./lever_arm_len; % Error in thrust

thrust_mean_EKF(kk) = norm(mean_torque_EKF)/lever_arm_len; % Find thrust

disp("EKF Measured thrust via omega: " + thrust_mean_EKF(kk).*1e6 + " +/- " + thrust_uncert_EKF.*1e6 + " uN")
%}
%% With Noise UKF Thrust diff Rotation
Tbb_UKF = zeros([3 length(t1)-1]);
Tbb_UKF_norms = zeros([1 length(t1)]);
for i = 2:length(t1)
    omega_b_b_2 = w_b_ECI_ukfilt(i,:)';
    omega_b_b_1 = w_b_ECI_ukfilt(i-1,:)';
    domega = (omega_b_b_2-omega_b_b_1)./delta_t; % Calculate change in omega, weighted by timestep
    Tbb_UKF(:,i) = J*domega+Functions.vec2cross(omega_b_b_1)*J*omega_b_b_1;
    
    quat_b_ECI_2 = quat_b_ECI_ukfilt(i,:)';
    C_b_ECI_2 = Functions.quat2rot(quat_b_ECI_2);
    Tbb_UKF(:,i) = C_b_ECI_2*Tbb_UKF(:,i);
    
    
    Tbb_UKF_norms(i) = norm(Tbb_UKF(:,i));

end
strIndex = 51; % at four seconds sim time

Tbb_UKF_calc = Tbb_UKF(:,strIndex:end); % Truncate torque values to approximately where convergence occurs
mean_torque_UKF = mean(Tbb_UKF_calc,2); % Find mean torque
mean_torque_uncert_UKF = std(Tbb_UKF_calc,0,2)./sqrt(length(t1)-1); % Calculate uncertainty in each component
norm_uncert_UKF = norm(mean_torque_uncert_UKF); % Error in torque vector norm is never greater than the norm of the errors
thrust_uncert_UKF = norm_uncert_UKF./lever_arm_len; % Error in thrust

thrust_mean_UKF(kk) = norm(mean_torque_UKF)/lever_arm_len; % Find thrust

disp("UKF Measured thrust via omega: " + thrust_mean_UKF(kk).*1e6 + " +/- " + thrust_uncert_UKF.*1e6 + " uN")
end

%% Plots

%%%%%%%
% EKF %
%%%%%%%
time_indx = zeros(size(t1));
for i = 1:length(t1)
    time_indx(i) = find(t1(i)==t);
end

err_states_EKF = w_b_ECI(time_indx,:)-w_b_ECI_ekfilt; % Find error in each estimated state
err_norm_EKF = vecnorm(err_states_EKF,2,2); % Perform the vector Euclidean norm on each row of the matrix

%%%%%%%
% UKF %
%%%%%%%

err_states_UKF = w_b_ECI(time_indx,:)-w_b_ECI_ukfilt; % Find error in each estimated state
err_norm_UKF = vecnorm(err_states_UKF,2,2); % Perform the vector Euclidean norm on each row of the matrix

%% Plot orbit
r = outKF.r.Data;

control_thrust_EKF = 23.284; % baseline measurement in micronewtons
control_thrust_UKF = 23.285; % baseline measurement in micronewtons

% Calculate standard deviation and range of values based on initial
% conditions variations
std_EKF = std(thrust_mean_EKF);
max_EKF = max(thrust_mean_EKF);
mean_EKF = mean(thrust_mean_EKF);
std_UKF = std(thrust_mean_UKF);
max_UKF = max(thrust_mean_UKF);
mean_UKF = mean(thrust_mean_UKF);

min_EKF = min(thrust_mean_EKF);
min_UKF = min(thrust_mean_UKF);

disp("EKF Thrust Ranges")
disp("Min: " + min_EKF.*1e6 + " uN")
disp("Max: " + max_EKF.*1e6 + " uN")
disp("Mean: " + mean_EKF.*1e6 + " uN")
disp("Std: " + std_EKF.*1e6 + " uN")

disp("UKF Thrust Ranges")
disp("Min: " + min_UKF.*1e6 + " uN")
disp("Max: " + max_UKF.*1e6 + " uN")
disp("Mean: " + mean_UKF.*1e6 + " uN")
disp("Std: " + std_UKF.*1e6 + " uN")

save('Spacecraft_outKF_Full_Sim_Vector_Uncert2')