%% Sub-Microthruster Fusion Technique Simulation
% Jonathan Hood, 1/8/2022
% Runs full simulation but performs Kalman Filter validation on both
% Kalman Filters.
% IMPORTANT NOTE: Uncomment magnetometer block in SIMULINK to verify validation OR
% modify variable dev_en to [1 0 1] to disable the magnetometer measurement
% analysis.
% Uncomment section in line 300 to generate validation data for analysis.

%% Spacecraft Geometry

clear
clc
close all

% CubeSat rectangular 6U design with thruster at one end perpendicular to center of
% mass of spacecraft. 
% Assume one thou for each directional tolerance
mech_tol = 0;
% Set spacecraft tolerancing values
sc_tolft = 1/1000;
sc_tol = 0.0254*sc_tolft;% spacecraft +/- tolerance in m

units = 6;               % number of cubes in configuration
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
thrust_loc_tol = 1/1000; % Set thruster location tolerance in in
thrust_tol = 0.0254*thrust_loc_tol*mech_tol;% thruster loc +/- tolerance in m

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
t1 = datetime(timeString1,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create MATLAB datetime object in PST

% Change t1 to UTC for math
JD_start = juliandate(t1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
y_sun = -cross(r,v)/norm(cross(r,v));              % y basis vector LVLH
x = cross(y_sun,z);                                % x basis vector LVLH

C_LVLH_ECI0 = [x;y_sun;z];                             % get ECI to LVLH rotation matrix 
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
drag_unmodeled = 1e-26.*drag_unmodeled_unnormalized./norm(drag_unmodeled_unnormalized).*dist_on;

% From sim, know SRP is essentially negligible (1e-24)
srp_unmodeled_unnormalized = randn([3 1]);
srp_unmodeled = 1e-24.*srp_unmodeled_unnormalized./norm(srp_unmodeled_unnormalized).*dist_on;

% From sim, know Mag is  nonnegligible (1e-5); must model to at least 1e-8
mag_field_unmodeled_unnormalized = randn([3 1]);
mag_field_unmodeled = 1e-9.*mag_field_unmodeled_unnormalized./norm(mag_field_unmodeled_unnormalized).*dist_on;

% From sim, know GG is nonnegligible (1e-7); must model to at least 1e-8,
% adjust as needed
earth_and_oblation_mag = 1e-9.*dist_on;
tidal_mag = 1e-9.*dist_on;
sun_and_moon_mag = 1e-14.*dist_on;

GG_mag = earth_and_oblation_mag+tidal_mag+sun_and_moon_mag;
gg_unmodeled_unnormalized = randn([3 1]);
grav_grad_unmodeled = GG_mag.*gg_unmodeled_unnormalized./norm(gg_unmodeled_unnormalized);

% General relativity
lens_thirring_mag = 1e-18.*dist_on;
schwarzchild_mag = 1e-15.*dist_on;
de_sitter_mag = 1e-18.*dist_on;

GR_mag = lens_thirring_mag+schwarzchild_mag+de_sitter_mag;
gr_unmodeled_unnormalized = randn([3 1]);
gr_eff_unmodeled = GG_mag.*gr_unmodeled_unnormalized./norm(gr_unmodeled_unnormalized);

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
sigma_g = 0.01; % gyroscope measurement noise
var_g = sigma_g^2; % gyroscope measurement noise variance

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

dev_en = [1 1 1]; % Determines which attitude determination devices will be considered in the sim run
                  % Use 0 to disable, 1 to enable

% Turn thrust on/off
thrust_on = 1;
% Set measurement timesteps
sun_meas_frequ = 1; % sun sensor measurement frequency in Hz
mag_meas_frequ = 1; % magnetometer measurement frequ in Hz
gyro_meas_frequ = 1; % gyroscope measurement frequ in Hz
g_Ts = 1/gyro_meas_frequ; % gyroscope sample time for gyro functions
ts = g_Ts;

meas_frequ_arr = [sun_meas_frequ mag_meas_frequ gyro_meas_frequ];
% Set integration timestep; must be integer divisor of measurement
% frequencies
max_frequ = max(meas_frequ_arr);
delta_t = 1/max_frequ; % time step between measurements & prediction time step
% Set simulation time
simtime = 20;       % desired simulation runtime in sec
% Set initial state guess
w_hat0 = [0;0;0]; % initial angular velocity
quat_hat0 = [0;0;0;1]; % initial quaternion
x_hat0 = [w_hat0; quat_hat0]; % initial state

%% Relative tolerance adjusted from default 1e-3 to 1e-12 to improve
% resolution
% theor_torque = theor_torque.*thrust_on;
% load_system("+Sims/attempt_2_Full_Microthruster_Simulation"); %load simulation
% set_param('attempt_2_Full_Microthruster_Simulation', 'StopTime', 'simtime') %set sim to stop at desired runtime
% outKF = sim('attempt_2_Full_Microthruster_Simulation');      % Run Simulink model
% save('outKF_KF_Validation_Data','outKF')
%% Data Analysis

load('outKF_KF_Validation_Data.mat')

% Extract continuous and discrete timedata
t = outKF.tout;
t1= outKF.x_hat_UKF.Time;

% Plot/Compare Measured Observation vs. Estimated Observation

% Extract continuous "real" measurements and reshape to 2D arrays
y_sun = outKF.y_sun.Data;
y_sun = reshape(y_sun,[3 length(y_sun)]);
y_mag = outKF.y_mag.Data;
y_mag = reshape(y_mag,[3 length(y_mag)]);
y_gyro = outKF.omega_gyro.Data';
% y_gyro = reshape(y_gyro,[3 length(y_gyro)]);

%%%%%%%%%%%%%%
%% Unscented %
%%%%%%%%%%%%%%

% Extract measurement data from EKF or UKF and reshape into 2D arrays
y_UKF_sun = outKF.y_UKF_sun.Data;
y_UKF_mag = outKF.y_UKF_mag.Data;
y_UKF_gyro = outKF.y_UKF_gyro.Data';

y_UKF_sun = reshape(y_UKF_sun,[3 length(y_UKF_sun)]);
y_UKF_mag = reshape(y_UKF_mag,[3 length(y_UKF_mag)]);
% y_UKF_gyro = reshape(y_UKF_gyro,[3 length(y_UKF_gyro)]);

UKF_sun_t = outKF.y_UKF_sun.Time;
UKF_sun_step = nnz(UKF_sun_t==0);
UKF_mag_t = outKF.y_UKF_mag.Time;
UKF_mag_step = nnz(UKF_mag_t==0);
UKF_gyro_t = outKF.y_UKF_gyro.Time;
UKF_gyro_step = nnz(UKF_gyro_t==0);


% Iterate through discrete timesteps and find indices that align with the
% continuous data
time_indx_sun = ones([1 length(UKF_sun_t)/UKF_sun_step]);
time_indx_mag = ones([1 length(UKF_mag_t)/UKF_mag_step]);
time_indx_gyro = ones([1 length(UKF_gyro_t)/UKF_gyro_step]);

count = 1;
for i = 1:UKF_sun_step:length(UKF_sun_t)
    indxs = find(UKF_sun_t(i)<=t);
    time_indx_sun(count) = indxs(1);
    count = count + 1;
end
count = 1;
for i = 1:UKF_mag_step:length(UKF_mag_t)
    indxs = find(UKF_mag_t(i)<=t);
    time_indx_mag(count) = indxs(1);
    count = count + 1;
end
count = 1;
for i = 1:UKF_gyro_step:length(UKF_gyro_t)
    indxs = find(UKF_gyro_t(i)<=t);
    time_indx_gyro(count) = indxs(1);
    count = count + 1;
end

% Extract measurements from continuous timesteps corresponding to discrete steps 
y_sun = y_sun(:,time_indx_sun);
y_mag = y_mag(:,time_indx_mag);
y_gyro = y_gyro(:,time_indx_gyro);

% Unsure on the procedure, but Simulink UKF block calls measurements
% functions a number of times affected by measurement frequency and UKF vs
% EKF procedure. Extract set elements
y_UKF_sun_trunc = y_UKF_sun(:,1:UKF_sun_step:end);
y_UKF_mag_trunc = y_UKF_mag(:,1:UKF_mag_step:end);
y_UKF_gyro_trunc = y_UKF_gyro(:,1:UKF_gyro_step:end);

% Calculate innovation for each measurement at each discrete timestep
Ik_sun = y_sun-y_UKF_sun_trunc;
Ik_mag = y_mag-y_UKF_mag_trunc;
Ik_gyro = y_gyro-y_UKF_gyro_trunc;

% Calculate autocorrelation
N_sun = length(time_indx_sun)-1;
N_mag = length(time_indx_mag)-1;
N_gyro = length(time_indx_gyro)-1;

rn_sun = zeros([3 (N_sun*2)+1]); % preallocate autocorrelation array
rn_mag = zeros([3 (N_mag*2)+1]); % preallocate autocorrelation array
rn_gyro = zeros([3 (N_gyro*2)+1]); % preallocate autocorrelation array
conf_int_sun = 1.96/sqrt(N_sun); % 95% confident it's only white noise, use 1.96/sqrt(N)
conf_int_mag = 1.96/sqrt(N_mag); % 95% confident it's only white noise, use 1.96/sqrt(N)
conf_int_gyro = 1.96/sqrt(N_gyro); % 95% confident it's only white noise, use 1.96/sqrt(N)

[rn_sun(1,:),lags_sun] = xcorr(Ik_sun(1,:),'normalized');
[rn_sun(2,:),lags_sun] = xcorr(Ik_sun(2,:),'normalized');
[rn_sun(3,:),lags_sun] = xcorr(Ik_sun(3,:),'normalized');

logi_sun = abs(rn_sun(:,N_sun+2:end))>=conf_int_sun;
sun_perc_out = nnz(logi_sun)/(N_sun*3)*100;

% Calculate autocorrelation of innovation for magnetometer
[rn_mag(1,:),lags_mag] = xcorr(Ik_mag(1,:),'normalized');
[rn_mag(2,:),lags_mag] = xcorr(Ik_mag(2,:),'normalized');
[rn_mag(3,:),lags_mag] = xcorr(Ik_mag(3,:),'normalized');

% Find all autocorrelations beyond 95% confidence interval and display
% percent outside; expect 5% or less
logi_mag = abs(rn_mag(:,N_mag+2:end))>=conf_int_mag;
mag_perc_out = nnz(logi_mag)/(N_mag*3)*100;


[rn_gyro(1,:),lags_gyro] = xcorr(Ik_gyro(1,:),'normalized');
[rn_gyro(2,:),lags_gyro] = xcorr(Ik_gyro(2,:),'normalized');
[rn_gyro(3,:),lags_gyro] = xcorr(Ik_gyro(3,:),'normalized');

logi_gyro = abs(rn_gyro(:,N_gyro+2:end))>=conf_int_gyro;
gyro_perc_out = nnz(logi_gyro)/(N_gyro*3)*100;


n = 1:N_sun;
figure
plot(n,rn_sun(1,N_sun+2:end),n,rn_sun(2,N_sun+2:end),n,rn_sun(3,N_sun+2:end))
hold on
grid on
xlabel('n-Step')
ylabel('Sun Autocorrelation')
yline(conf_int_sun,'--k')
yline(-conf_int_sun,'--k')
title('Unscented Kalman Filter Sun Sensor Autocorrelation Plot')
legend('Sun 1R_n','Sun 2R_n','Sun 3R_n',...
    '95% Conf. Int.')
exportgraphics(gcf,'UKF_sun_autocorrelation.png')

n = 1:N_mag;
figure
plot(n,rn_mag(1,N_mag+2:end),n,rn_mag(2,N_mag+2:end),n,rn_mag(3,N_mag+2:end))
hold on
grid on
xlabel('n-Step')
ylabel('Mag Autocorrelation')
yline(conf_int_sun,'--k')
yline(-conf_int_sun,'--k')
title('Unscented Kalman Filter Magnetometer Autocorrelation Plot')
legend('Mag 1R_n','Mag 2R_n','Mag 3R_n',...
    '95% Conf. Int.')
exportgraphics(gcf,'UKF_mag_autocorrelation.png')

n = 1:N_gyro;
figure
plot(n,rn_gyro(1,N_gyro+2:end),n,rn_gyro(2,N_gyro+2:end),n,rn_gyro(3,N_gyro+2:end))
hold on
grid on
xlabel('n-Step')
ylabel('Gyro Autocorrelation')
yline(conf_int_gyro,'--k')
yline(-conf_int_gyro,'--k')
title('Unscented Kalman Filter Gyroscope Autocorrelation Plot')
legend('Gyro 1R_n','Gyro 2R_n','Gyro 3R_n',...
    '95% Conf. Int.')
exportgraphics(gcf,'UKF_gyro_autocorrelation.png')


disp("=======")
disp("= UKF =")
disp("=======")
fprintf('\nMax mag autocorrelation UKF: %e\n',max(max(abs(rn_mag(2,N_mag+2:end)))))
fprintf('UFK magnetometer 95%% confidence interval: %f\n',conf_int_mag)
fprintf('UFK magnetometer %f%% outside confidence\n',mag_perc_out)

fprintf('\nMax sun autocorrelation UKF: %e\n',max(max(abs(rn_sun(1,N_sun+2:end)))))
fprintf('UFK sun sensor 95%% confidence interval: %f\n',conf_int_sun)
fprintf('UFK sun sensor %f%% outside confidence\n',sun_perc_out)

fprintf('\nMax gyro autocorrelation UKF: %e\n',max(max(abs(rn_gyro(1,N_gyro+2:end)))))
fprintf('UFK gyro sensor 95%% confidence interval: %f\n',conf_int_gyro)
fprintf('UFK gyro sensor %f%% outside confidence\n',gyro_perc_out)

%%%%%%%%%%%%%%
%% Extended  %
%%%%%%%%%%%%%%

% Extract continuous "real" measurements and reshape to 2D arrays
y_sun = outKF.y_sun.Data;
y_sun = reshape(y_sun,[3 length(y_sun)]);
y_mag = outKF.y_mag.Data;
y_mag = reshape(y_mag,[3 length(y_mag)]);
y_gyro = outKF.omega_gyro.Data';
% y_gyro = reshape(y_gyro,[3 length(y_gyro)]);

% Extract measurement data from EKF or UKF and reshape into 2D arrays
y_EKF_sun = outKF.y_EKF_sun.Data;
y_EKF_mag = outKF.y_EKF_mag.Data;
y_EKF_gyro = outKF.y_EKF_gyro.Data';

y_EKF_sun = reshape(y_EKF_sun,[3 length(y_EKF_sun)]);
y_EKF_mag = reshape(y_EKF_mag,[3 length(y_EKF_mag)]);
% y_EKF_gyro = reshape(y_EKF_gyro,[3 length(y_EKF_gyro)]);

EKF_sun_t = outKF.y_EKF_sun.Time;
EKF_sun_step = nnz(EKF_sun_t==0);
EKF_mag_t = outKF.y_EKF_mag.Time;
EKF_mag_step = nnz(EKF_mag_t==0);
EKF_gyro_t = outKF.y_EKF_gyro.Time;
EKF_gyro_step = nnz(EKF_gyro_t==0);

% Iterate through discrete timesteps and find indices that align with the
% continuous data
time_indx_sun = ones([1 length(EKF_sun_t)/EKF_sun_step]);
time_indx_mag = ones([1 length(EKF_mag_t)/EKF_mag_step]);
time_indx_gyro = ones([1 length(EKF_gyro_t)/EKF_gyro_step]);
count = 1;
for i = 1:EKF_sun_step:length(EKF_sun_t)
    indxs = find(EKF_sun_t(i)<=t);
    time_indx_sun(count) = indxs(1);
    count = count + 1;
end
count = 1;
for i = 1:EKF_mag_step:length(EKF_mag_t)
    indxs = find(EKF_mag_t(i)<=t);
    time_indx_mag(count) = indxs(1);
    count = count + 1;
end
count = 1;
for i = 1:EKF_gyro_step:length(EKF_gyro_t)
    indxs = find(EKF_gyro_t(i)<=t);
    time_indx_gyro(count) = indxs(1);
    count = count + 1;
end

% Extract measurements from continuous timesteps corresponding to discrete steps 
y_sun = y_sun(:,time_indx_sun);
y_mag = y_mag(:,time_indx_mag);
y_gyro = y_gyro(:,time_indx_gyro);

% Unsure on the procedure, but Simulink EKF block calls measurements
% functions a number of times affected by measurement frequency and EKF vs
% EKF procedure. Extract set elements
y_EKF_sun_trunc = y_EKF_sun(:,1:EKF_sun_step:end);
y_EKF_mag_trunc = y_EKF_mag(:,1:EKF_mag_step:end);
y_EKF_gyro_trunc = y_EKF_gyro(:,1:EKF_gyro_step:end);

% Calculate innovation for each measurement at each discrete timestep
Ik_sun = y_sun-y_EKF_sun_trunc;
Ik_mag = y_mag-y_EKF_mag_trunc;
Ik_gyro = y_gyro-y_EKF_gyro_trunc;

% Calculate autocorrelation
N_sun = length(time_indx_sun)-1;
N_mag = length(time_indx_mag)-1;
N_gyro = length(time_indx_gyro)-1;

rn_sun = zeros([3 (N_sun*2)+1]); % preallocate autocorrelation array
rn_mag = zeros([3 (N_mag*2)+1]); % preallocate autocorrelation array
rn_gyro = zeros([3 (N_gyro*2)+1]); % preallocate autocorrelation array
conf_int_sun = 1.96/sqrt(N_sun); % 95% confident it's only white noise, use 1.96/sqrt(N)
conf_int_mag = 1.96/sqrt(N_mag); % 95% confident it's only white noise, use 1.96/sqrt(N)
conf_int_gyro = 1.96/sqrt(N_gyro); % 95% confident it's only white noise, use 1.96/sqrt(N)

[rn_sun(1,:),lags_sun] = xcorr(Ik_sun(1,:),'normalized');
[rn_sun(2,:),lags_sun] = xcorr(Ik_sun(2,:),'normalized');
[rn_sun(3,:),lags_sun] = xcorr(Ik_sun(3,:),'normalized');

logi_sun = abs(rn_sun(:,N_sun+2:end))>=conf_int_sun;
sun_perc_out = nnz(logi_sun)/(N_sun*3)*100;

% Calculate autocorrelation of innovation for magnetometer
[rn_mag(1,:),lags_mag] = xcorr(Ik_mag(1,:),'normalized');
[rn_mag(2,:),lags_mag] = xcorr(Ik_mag(2,:),'normalized');
[rn_mag(3,:),lags_mag] = xcorr(Ik_mag(3,:),'normalized');

% Find all autocorrelations beyond 95% confidence interval and display
% percent outside; expect 5% or less
logi_mag = abs(rn_mag(:,N_mag+2:end))>=conf_int_mag;
mag_perc_out = nnz(logi_mag)/(N_mag*3)*100;


% Calculate autocorrelation of innovation for gyroscope
[rn_gyro(1,:),lags_gyro] = xcorr(Ik_gyro(1,:),'normalized');
[rn_gyro(2,:),lags_gyro] = xcorr(Ik_gyro(2,:),'normalized');
[rn_gyro(3,:),lags_gyro] = xcorr(Ik_gyro(3,:),'normalized');

% Find all autocorrelations beyond 95% confidence interval and display
% percent outside; expect 5% or less
logi_gyro = abs(rn_gyro(:,N_gyro+2:end))>=conf_int_gyro;
gyro_perc_out = nnz(logi_gyro)/(N_gyro*3)*100;


n = 1:N_sun;
figure
plot(n,rn_sun(1,N_sun+2:end),n,rn_sun(2,N_sun+2:end),n,rn_sun(3,N_sun+2:end))
hold on
grid on
xlabel('n-Step')
ylabel('Sun Autocorrelation')
yline(conf_int_sun,'--k')
yline(-conf_int_sun,'--k')
title('Extended Kalman Filter Sun Sensor Autocorrelation Plot')
legend('Sun 1R_n','Sun 2R_n','Sun 3R_n',...
    '95% Conf. Int.')
exportgraphics(gcf,'EKF_sun_autocorrelation.png')

n = 1:N_mag;
figure
plot(n,rn_mag(1,N_mag+2:end),n,rn_mag(2,N_mag+2:end),n,rn_mag(3,N_mag+2:end))
hold on
grid on
xlabel('n-Step')
ylabel('Mag Autocorrelation')
yline(conf_int_sun,'--k')
yline(-conf_int_sun,'--k')
title('Extended Kalman Filter Magnetometer Autocorrelation Plot')
legend('Mag 1R_n','Mag 2R_n','Mag 3R_n',...
    '95% Conf. Int.')
exportgraphics(gcf,'EKF_mag_autocorrelation.png')

n = 1:N_gyro;
figure
plot(n,rn_gyro(1,N_gyro+2:end),n,rn_gyro(2,N_gyro+2:end),n,rn_gyro(3,N_gyro+2:end))
hold on
grid on
xlabel('n-Step')
ylabel('Gyro Autocorrelation')
yline(conf_int_sun,'--k')
yline(-conf_int_sun,'--k')
title('Extended Kalman Filter Gyroscope Autocorrelation Plot')
legend('Gyro 1R_n','Gyro 2R_n','Gyro 3R_n',...
    '95% Conf. Int.')
exportgraphics(gcf,'EKF_gyro_autocorrelation.png')


disp("=======")
disp("= EKF =")
disp("=======")
fprintf('\nMax mag autocorrelation EKF: %e\n',max(max(abs(rn_mag(2,N_mag+2:end)))))
fprintf('EFK magnetometer 95%% confidence interval: %f\n',conf_int_mag)
fprintf('EFK magnetometer %f%% outside confidence\n',mag_perc_out)

fprintf('\nMax sun autocorrelation EKF: %e\n',max(max(abs(rn_sun(1,N_sun+2:end)))))
fprintf('EFK sun sensor 95%% confidence interval: %f\n',conf_int_sun)
fprintf('EFK sun sensor %f%% outside confidence\n',sun_perc_out)

fprintf('\nMax gyro autocorrelation EKF: %e\n',max(max(abs(rn_gyro(1,N_gyro+2:end)))))
fprintf('EFK gyro sensor 95%% confidence interval: %f\n',conf_int_gyro)
fprintf('EFK gyro sensor %f%% outside confidence\n',gyro_perc_out)
