%% Sub-Microthruster Fusion Technique Simulation: Orbit Adjustment
% Jonathan Hood, 5/3/2022
% Simulate an on-orbit sub-microthruster evaluation via the proposed Fusion
% Technique, including orbital dynamics and present perturbations. Run
% through a range of varying inclinations and orbital altitudes to quantify
% effect on measured thrust.
%% Spacecraft Geometry
% To be modified and iterated to meet performance requirements (10^-7 or >
% thrust resolution)

clear
clc
close all

% Assume one thou for each directional tolerance
mech_tol = 0;
% Set spacecraft tolerancing values
sc_tol = 1/1000;% spacecraft +/- tolerance in m
orb_range = [300 400 550 2000 5000 10000 20000 30000 40000];
inc_range = [0 20 40 60 80 100 120 140 160 180];

for kk = 1:length(orb_range)
    for jj = 1:length(inc_range)
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
thrust_tol = 1/1000*mech_tol;% thruster loc +/- tolerance in m

lever_arm_len = L/2; % lever arm length in m
theor_thrust = 1e-7; % predicted value of desired thrust to be measured
thruster_loc = [L/2+randn(1)*thrust_tol;thrust_tol*randn(1);thrust_tol*randn(1)]; % location of thruster at far edge of x-axis

az_tol = randn(1)*0.001*mech_tol; % Error in thruster force vector in the xy plane
dec_tol = randn(1)*0.001*mech_tol; % Error in thruster force vector about the z-plane
theta = 0; % optimal azimuthal theta is aligned with the y-axis
gamma = 0; % optimal declination gamma is 0, or in the xy plane
thrust_vec_unit = [sind(theta+az_tol); cosd(az_tol+theta); sind(dec_tol+gamma)];
thrust_vec_unit = thrust_vec_unit./norm(thrust_vec_unit);
thrust_vec = theor_thrust.*thrust_vec_unit; % Constant thrust vector in the positive y-dir
theor_torque = cross(thruster_loc,thrust_vec); % torque vector is cross of thruster position and thrust vector

%% Orbital Motion Parameters
% Assumed to start in a rideshare SSO orbit.

% Calculate position
alt = orb_range(kk);          % s/c altitude in km (from SpaceX average rideshare SSO altitude)
rE = 6378;          % radius of the Earth in km SPECIFIC
mu = 398600;        % std. grav param km^3/s^2
rmag = rE+alt;      % s/c position in km
inc = inc_range(jj);           % relatively ARBITRARY SSO inclination in deg
ecc = 0;            % assume circular orbit
h = sqrt(rmag*mu);  % relative angular momentum of circular
RA = 0;             % let RAAN be 0
w = 0;              % let argument of perigee be 0
TA = 0;             % let TA be 0

coe = [h, ecc, RA, inc, w, TA];         % Create COES vector
[r,v] = Functions.sv_from_coe(coe,mu);  % get initial r & v vector in km and km/s
COES = Functions.COES(r,v,mu);          % Extract more detailed COES
% AERO 557 HW-5 NCKF Unit Test of Orbital Motion Model
% Uncomment the below and save the position data for Orbital Motion
% validation. Compare 'r' results to 'r' results found running "HW5.m".
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = diag([27 17 25]); % s/c mass moment of inertia kg*m^2
simtime = 400; % sim in s

% Calculate position
alt = 450;          % s/c altitude in km
rE = 6378;          % radius of the Earth in km
mu = 398600;        % std. grav param
rmag = rE+alt;      % s/c position
inc = 87;           % inclination in deg
ecc = 0;            % circular orbit
h = sqrt(rmag*mu);  % relative angular momentum
RA = 0;             % let RAAN be 0
w = 0;              % let argument of perigee be 0
TA = 0;             % let TA be 0

coe = [h, ecc, RA, inc, w, TA]; % Create COES vecotr
[r,v] = Functions.sv_from_coe(coe,mu); % get r & v vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

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

% w_b_ECI0 = [0.0; 0; 0]; % initial angular velocity rad/s ARBITRARY from detumble
w_b_ECI0 = [0.05; -0.05; 0.05]; % initial angular velocity rad/s
quat_b_ECI0 = [sin(1/2)/sqrt(3)*ones([3 1]); cos(1/2)]; % initial quaternion ARBITRARY from detumble
init_state = [w_b_ECI0; quat_b_ECI0];

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

% % TEST initial covariance
% P0 = diag([0.1 0.1 0.1]);


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

% M_k = eye(size(R_k)); % identity matrix size of R_k

R = diag([sigma_m^2,sigma_m^2,sigma_m^2,sigma_s^2,sigma_s^2,sigma_s^2]); % initial measurement noise covariance
R_norm = norm(R);

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
tic
outKF = sim('attempt_2_Full_Microthruster_Simulation');      % Run Simulink model
time_out = toc;
% save('outKF_Full_Sim_no_mag')

% load('outKF_Full_Sim_control.mat')
fprintf('Simulation time %f s\n',time_out)

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

if length(sun_t) > length(gyro_t)
    t11 = sun_t;
else
    t11 = gyro_t;
end

UKF_step = nnz(t11==0);
time_indx = ones([1 length(t11)/UKF_step]);

count = 1;
for i = 1:UKF_step:length(t11)
    indxs = find(t11(i)<=t);
    time_indx(count) = indxs(1);
    count = count + 1;
end

t1 = t(time_indx);

%% No Noise True Thrust diff Rotation

Tbb_true = zeros([3 length(t)]);
Tbb_true_norms = zeros([1 length(t)]);
for i = 1:length(t)
    domega = domega_b_ECI(i,:)'; % change in angular velocity
    omega_b_b_1 = w_b_ECI(i,:)'; % angular velocity at that timestep
    Tbb_true(:,i) = J*domega+Functions.vec2cross(omega_b_b_1)*J*omega_b_b_1;
    quat_b_ECI_2 = quat_b_ECI(i,:)';
    C_b_ECI_2 = Functions.quat2rot(quat_b_ECI_2);
    Tbb_true(:,i) = C_b_ECI_2*Tbb_true(:,i);
Tbb_true_norms(i) = norm(Tbb_true(:,i));
end
mean_torque_true = mean(Tbb_true,2); % Find mean torque
mean_torque_uncert = std(Tbb_true,0,2)./sqrt(length(t)); % Calculate uncertainty in each component
norm_uncert = norm(mean_torque_uncert); % Error in torque vector norm is never greater than the norm of the errors
thrust_uncert = norm_uncert./lever_arm_len; % Error in thrust

thrust_mean = norm(mean_torque_true)/lever_arm_len; % Find thrust

% Check axis-specific initial angular velocity effects
disp("Mean true thrust: " + thrust_mean*1e6 + " +/- " + thrust_uncert*1e6 + " N")
figure
plot(t,Tbb_true_norms)
grid on
xlabel('Time (s)')
ylabel('Torque Magnitude (N)')
exportgraphics(gcf,'control_torque_mag.png')
Tbb_true_one = J*domega+Functions.vec2cross(omega_b_b_1)*J*omega_b_b_1;
thrust_meas = norm(Tbb_true_one)/lever_arm_len; % Since torque is now perpendicular to r, T = r*F
disp("Actual thrust: " + thrust_meas*1e6 + " uN")

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

figure
plot(t1(strIndex:end),Tbb_EKF_norms(strIndex:end))
grid on
xlabel('Time (s)')
ylabel('Torque Magnitude: EKF (N)')

Tbb_EKF_calc = Tbb_EKF(:,strIndex:end); % Truncate torque values to approximately where convergence occurs
mean_torque_EKF = mean(Tbb_EKF_calc,2); % Find mean torque
mean_torque_uncert_EKF = std(Tbb_EKF_calc,0,2)./sqrt(length(Tbb_EKF_calc)); % Calculate uncertainty in each component
norm_uncert_EKF = norm(mean_torque_uncert_EKF); % Error in torque vector norm is never greater than the norm of the errors
thrust_uncert_EKF = norm_uncert_EKF./lever_arm_len; % Error in thrust

thrust_mean_EKF(kk,jj) = norm(mean_torque_EKF)/lever_arm_len; % Find thrust

disp("EKF Measured thrust via omega: " + thrust_mean_EKF(kk,jj).*1e6 + " +/- " + thrust_uncert_EKF.*1e6 + " uN")
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

figure
plot(t1(strIndex:end),Tbb_UKF_norms(strIndex:end))
grid on
xlabel('Time (s)')
ylabel('Torque Magnitude: UKF (N)')

Tbb_UKF_calc = Tbb_UKF(:,strIndex:end); % Truncate torque values to approximately where convergence occurs
mean_torque_UKF = mean(Tbb_UKF_calc,2); % Find mean torque
mean_torque_uncert_UKF = std(Tbb_UKF_calc,0,2)./sqrt(length(t1)-1); % Calculate uncertainty in each component
norm_uncert_UKF = norm(mean_torque_uncert_UKF); % Error in torque vector norm is never greater than the norm of the errors
thrust_uncert_UKF = norm_uncert_UKF./lever_arm_len; % Error in thrust

thrust_mean_UKF(kk,jj) = norm(mean_torque_UKF)/lever_arm_len; % Find thrust

disp("UKF Measured thrust via omega: " + thrust_mean_UKF(kk,jj).*1e6 + " +/- " + thrust_uncert_UKF.*1e6 + " uN")
close all
    end
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
% Plot norm of error over time
f_EKF_err = figure;
plot(t1,err_norm_EKF)
hold on
grid on
xlabel('Time (s)')
ylabel('\omega Error Norm (rad/s)')
title('Angular Velocity: EKF Error')

% Plot angular velocity
f_EKF = figure;
subplot(1,2,1)
plot(t,w_b_ECI(:,1),t,w_b_ECI(:,2),t,w_b_ECI(:,3))
hold on
grid on
xlabel('Time (s)')
ylabel('\omega (rad/s)')
title('Angular Velocity: True vs. EKF')
plot(t1,w_b_ECI_ekfilt(:,1),'--',t1,w_b_ECI_ekfilt(:,2),'--',...
   t1,w_b_ECI_ekfilt(:,3),'--')
legend('\omega_1','\omega_2','\omega_3','\omega_{EKF1}','\omega_{EKF2}',...
    '\omega_{EKF3}')

% Plot quaternion components
subplot(1,2,2)
plot(t,quat_b_ECI(:,1),t,quat_b_ECI(:,2),t,quat_b_ECI(:,3),...
    t,quat_b_ECI(:,4))
hold on
grid on
xlabel('Time (s)')
ylabel('Magnitude')
title('Quaternion: True vs. EKF')
plot(t1,quat_b_ECI_ekfilt(:,1),'--',t1,quat_b_ECI_ekfilt(:,2),'--',...
   t1,quat_b_ECI_ekfilt(:,3),'--',t1,quat_b_ECI_ekfilt(:,4),'--')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta','\epsilon_{EKF1}',...
    '\epsilon_{EKF2}','\epsilon_{EKF3}','\eta_{EKF}')

f_EKF.Units = 'normalized';
f_EKF.Position = [0 0.3 0.5 0.5];

exportgraphics(gcf,'ekf_filter_results.png')
%}
%%%%%%%
% UKF %
%%%%%%%

err_states_UKF = w_b_ECI(time_indx,:)-w_b_ECI_ukfilt; % Find error in each estimated state
err_norm_UKF = vecnorm(err_states_UKF,2,2); % Perform the vector Euclidean norm on each row of the matrix
% Plot norm of error over time
f_UKF_err = figure;
plot(t1,err_norm_UKF)
hold on
grid on
xlabel('Time (s)')
ylabel('\omega Error Norm (rad/s)')
title('Angular Velocity: UKF Error')

f_UKF = figure;
subplot(1,2,1)
plot(t,w_b_ECI(:,1),t,w_b_ECI(:,2),t,w_b_ECI(:,3))
hold on
grid on
xlabel('Time (s)')
ylabel('\omega (rad/s)')
title('Angular Velocity: True vs. UKF')
plot(t1,w_b_ECI_ukfilt(:,1),'--',t1,w_b_ECI_ukfilt(:,2),'--',...
   t1,w_b_ECI_ukfilt(:,3),'--')
legend('\omega_1','\omega_2','\omega_3','\omega_{UKF1}','\omega_{UKF2}',...
    '\omega_{UKF3}')


% 1.96/sqrt(N) for 95 conf int; 
% See how variance/uncertainty varies with input sampling rates,
% measurement noise, sensor type, simulation time (ideally small), geometry
% Remeber to limit resolution of attitude sensors
% 

% Plot quaternion components
subplot(1,2,2)
plot(t,quat_b_ECI(:,1),t,quat_b_ECI(:,2),t,quat_b_ECI(:,3),...
    t,quat_b_ECI(:,4))
hold on
grid on
xlabel('Time (s)')
ylabel('Magnitude')
title('Quaternion: True vs. UKF')
plot(t1,quat_b_ECI_ukfilt(:,1),'--',t1,quat_b_ECI_ukfilt(:,2),'--',...
   t1,quat_b_ECI_ukfilt(:,3),'--',t1,quat_b_ECI_ukfilt(:,4),'--')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta','\epsilon_{UKF1}',...
    '\epsilon_{UKF2}','\epsilon_{UKF3}','\eta_{UKF}')

f_UKF.Units = 'normalized';
f_UKF.Position = [0.5 0.3 0.5 0.5];

exportgraphics(gcf,'ukf_filter_results.png')


%% Plot orbit
r = outKF.r.Data;

% Save r data
simtime_str = sprintf('%.2f',simtime/60);
save('r_micro','r')
Functions.earth
rfig = gcf;
hold on
scatter3(r(:,1),r(:,2),r(:,3),'filled')
xlabel('X (km)')
scatter3(r(1,1),r(1,2),r(1,3),50,"green",'filled')
scatter3(r(end,1),r(end,2),r(end,3),50,"red",'filled')
ylabel('Y (km)')
zlabel('Z (km)')
title("Orbital Motion over Simtime = " + simtime_str + " min in ECI Frame")
legend('Orbit','Start','End')
view(70,18)
rfig.Units = 'normalized';
rfig.Position = [0.4 0.3 0.3 0.6];
exportgraphics(gcf,'orbital_motion.png')
control_thrust_EKF = 23.284; % baseline measurement in micronewtons
control_thrust_UKF = 23.285; % baseline measurement in micronewtons

% Calculate and display the reduction in measured thrust based on orbital
% parameters
for kk = 1:length(orb_range)
    for jj = 1:length(inc_range)
    disp("EKF Measured thrust via omega: " + thrust_mean_EKF(kk,jj).*1e6 + "uN")
    disp("UKF Measured thrust via omega: " + thrust_mean_UKF(kk,jj).*1e6 + "uN")
    EKF_perc_reduc = 1-(thrust_mean_EKF(kk,jj)*1e6/control_thrust_EKF);
    UKF_perc_reduc = 1-(thrust_mean_UKF(kk,jj)*1e6/control_thrust_UKF);
    if thrust_mean_EKF(kk,jj)*1e6 < control_thrust_EKF
        fprintf('Alt %d, Inc %d, EKF Results reduced by %f %%\n',orb_range(kk),inc_range(jj),EKF_perc_reduc*100)
    end
    if thrust_mean_UKF(kk,jj)*1e6 < control_thrust_UKF
        fprintf('\nAlt %d, Inc %d, UKF results reduced by %f %%\n\n',orb_range(kk),inc_range(jj),UKF_perc_reduc*100)
    end
    end
end

save('Orbit_outKF_Full_Sim')