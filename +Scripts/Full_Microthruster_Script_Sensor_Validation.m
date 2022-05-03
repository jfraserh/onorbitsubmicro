%% Sub-Microthruster Fusion Technique Simulation: Sensor Validation and Perturbation Analysis
% Jonathan Hood, 5/3/2022
% Runs full simulation but performs validation on ADDs; also performs
% perturbation analysis.
% Modification to include nested for loops also incorporates perturbation
% analysis at the included orb_range and inc_range values.

%% Spacecraft Geometry
% To be modified and iterated to meet performance requirements (10^-7 or >
% thrust resolution)

clear
clc
close all

% CubeSat rectangular 6U design with thruster at one end perpendicular to center of
% mass of spacecraft. 
% Assume one thou for each directional tolerance
mech_tol = 0;
% Set spacecraft tolerancing values
sc_tol = 1/1000;% spacecraft +/- tolerance in m
orb_range = [300 400 550 2000 5000 10000 20000 30000 40000];
inc_range = [0 20 40 60 80 100 120 140 160 180];

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
% AERO 557 HW-5 NCKF Unit Test of Orbital Motion Model

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
z1 = -r/norm(r);                                % z basis vector LVLH
y1 = -cross(r,v)/norm(cross(r,v));              % y basis vector LVLH
x1 = cross(y1,z1);                                % x basis vector LVLH

C_LVLH_ECI0 = [x1;y1;z1];                             % get ECI to LVLH rotation matrix 
quat_LVLH_ECI0 = Functions.rot2quat(C_LVLH_ECI0);   % get ECI to LVLH quaternion
euler_LVLH_ECI0 = Functions.rot2euler(C_LVLH_ECI0); % get ECI to LVLH euler angles

% Initial b_LVLH
C_b_LVLH0 = C_b_ECI0*C_LVLH_ECI0';                     % get b to LVLH rotation matrix 
quat_b_LVLH0 = Functions.rot2quat(C_b_LVLH0);         % get b to LVLH quaternion
euler_b_LVLH0 = Functions.rot2euler(C_b_LVLH0);       % get b to LVLH euler angles

%% Sun Torque Setup

% Extract s/c geometry parameters
x = geom(1); % x-axis
y = geom(3); % y-axis
z = geom(2); % z-axis

%%%%%% ATMO DRAG %%%%%%%

% Instantiate CoP vectors
CoM_vec = [0;0;0]; % From mass_props inertia calculations; assume CoM is at center of body frame

%%%%%%%%%BUS COP VECTORS%%%%%%%%%%%%%%%%
BUS_neg_z = [0;0;-z/2];
CP_BUS_neg_z = CoM_vec - BUS_neg_z;

BUS_pos_z = -BUS_neg_z;
CP_BUS_pos_z = CoM_vec - BUS_pos_z;

BUS_neg_y = [0;-y/2;0];
CP_BUS_neg_y = CoM_vec - BUS_neg_y;

BUS_pos_y = -BUS_neg_y;
CP_BUS_pos_y = CoM_vec - BUS_pos_y;

BUS_neg_x = [-x/2;0;0];
CP_BUS_neg_x = CoM_vec - BUS_neg_x;

BUS_pos_x = -BUS_neg_x;
CP_BUS_pos_x = CoM_vec - BUS_pos_x;

% Find area of each exposed surface
A_SP_y_side = x*z;
A_SP_x_side = y*z;
A_SP_z_side = x*y;

% Create normal directional vectors
x_dir = [1;0;0];
y_dir = [0;1;0];
z_dir = [0;0;1];
dir_vector = {x_dir,-x_dir,y_dir,-y_dir,z_dir,-z_dir};

% Create arrays of areas and CoM-->CoP vectors for each surface for each
% normal direction (x, -x, y, -y, z, -z)
pos_x_COP_vectors = [CP_BUS_pos_x];
pos_x_areas       = [A_SP_x_side];

neg_x_areas        = pos_x_areas;
neg_x_COP_vectors = [CP_BUS_neg_x];

pos_y_COP_vectors = [CP_BUS_pos_y];
pos_y_areas        = [A_SP_y_side];

neg_y_areas        = pos_y_areas;
neg_y_COP_vectors = [CP_BUS_neg_y];

neg_z_areas        = [A_SP_z_side];
neg_z_COP_vectors = [CP_BUS_neg_z];

pos_z_areas        = [A_SP_z_side];
pos_z_COP_vectors = [CP_BUS_pos_z];

% Combine cell arrays of area and their corresponding vectors into parent
% cell arrays
CoP_vectors = {pos_x_COP_vectors,neg_x_COP_vectors,...
               pos_y_COP_vectors,neg_y_COP_vectors,...
               pos_z_COP_vectors,neg_z_COP_vectors};
areas = {pos_x_areas,neg_x_areas,pos_y_areas,neg_y_areas,pos_z_areas,neg_z_areas};

p_s = 4.5e-06; % solar const N/m^2
Tsrp = [0;0;0];                    % Initialzie atmo torque vector


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
% L_k = [delta_t*inv(J); zeros(3); zeros([1 3])];

% Sun and Magnetometer Measurement vectors
% Calculate initial sun measurement vector in ECI via Algorithm 4.3 in Hall
% 2003, Ch.4 Attitude Determination. 

jul_cent = (JD_start-2451545)/36525; % julian centuries since 2000
lambda_S = wrapTo360(280.4606184 + 36000.77005361*jul_cent); % mean longitude of sun in deg
M_S = wrapTo360(357.5277233 + 35999.05034*jul_cent); % mean anomaly of sun in deg
lambda_e = wrapTo180(lambda_S + 1.914666471*sind(M_S) + 0.019994643*sind(2*M_S)); % ecliptic long of sun in deg

eps = 23.439291-0.0130042*jul_cent;
sun_ECI0 = [cosd(lambda_e); cosd(eps)*sind(lambda_e); sind(eps)*sind(lambda_e)];
sun_hat_ECI0 = sun_ECI0./norm(sun_ECI0);

disp("==================")
disp("= Sun Validation =")
disp("==================")
JPL_vec = [1.098785289384590e8  -9.033968046078207e7  -3.916233990143497e7];
JPL_vec_unit = JPL_vec./norm(JPL_vec);
fprintf('\nJPL Ephemerides Vector ECI:  %f %f %f\n',JPL_vec_unit(1),JPL_vec_unit(2),JPL_vec_unit(3))
fprintf('Hall Algorithm Results ECI: %f %f %f\n',sun_hat_ECI0(1),sun_hat_ECI0(2),sun_hat_ECI0(3))

[y, m, d] = ymd(t1);
[h, mm, s] = hms(t1);
utc_t1 = [y m d h mm s];

disp("==================")
disp("= Mag Validation =")
disp("==================")

igrf_m_NED = igrfmagm(alt_g,del,lambda,decyear(y,m,d))*1e-9;
igrf_m_mag = norm(igrf_m_NED);
[igrf_m_ECEF(1),igrf_m_ECEF(2),igrf_m_ECEF(3)] = ned2ecefv(igrf_m_NED(1),igrf_m_NED(2),igrf_m_NED(3),del,lambda);
C_ECEF_ECI = dcmeci2ecef('IAU-2000/2006',utc_t1);
igrf_m_ECI = C_ECEF_ECI'*igrf_m_ECEF';

IGRF_NED = [23574.4 4015.1 2715.8]*1e-9;
% https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml?model=igrf#igrfwmm
% NOAA calculated values for given long, lat, alt
fprintf('\nOnline IGRF model vector in NED:  %e %e %e\n',IGRF_NED(1),IGRF_NED(2),IGRF_NED(3))
fprintf('MATLAB IGRF-13 Mag Field Vector in NED: %e %e %e\n',igrf_m_NED(1),igrf_m_NED(2),igrf_m_NED(3))

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
sun_meas_frequ = 1; % sun sensor measurement frequency in Hz
mag_meas_frequ = 1; % magnetometer measurement frequ in Hz
gyro_meas_frequ = 50; % gyroscope measurement frequ in Hz
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
%% Only run the sim when necessary
load_system("+Sims/attempt_2_Full_Microthruster_Simulation"); %load simulation
set_param('attempt_2_Full_Microthruster_Simulation', 'StopTime', 'simtime') %set sim to stop at desired runtime
outKF = sim('attempt_2_Full_Microthruster_Simulation');      % Run Simulink model
% save('outKF_Sensor_Validation','outKF');
%% Load previous sim run data for comparison
% % Extract continuous and discrete timedata

% load outKF_1U_no_mech.mat

t = outKF.tout;
t1= outKF.x_hat_UKF.Time;
y_UKF_sun = outKF.y_UKF_sun.Data;
y_UKF_mag = outKF.y_UKF_mag.Data;
y_UKF_gyro = outKF.y_UKF_gyro.Data;
w_b_ECI = outKF.w_b_ECI.Data;
domega_b_ECI = outKF.domega_b_ECI.Data; % Extract differential angular velocity for torque calculations
x_hat_EKF = outKF.x_hat_EKF.Data;
x_hat_UKF = outKF.x_hat_UKF.Data;
gyro_no_noise = outKF.omega_gyro_no_noise.Data;

% x_hat_k = outKF.x_hat_k.Data;
w_b_ECI_ekfilt = x_hat_EKF(:,1:3);
quat_b_ECI_ekfilt = x_hat_EKF(:,4:7);
w_b_ECI_ukfilt = x_hat_UKF(:,1:3);
quat_b_ECI_ukfilt = x_hat_UKF(:,4:7);
quat_b_ECI = outKF.quats_b_ECI.Data;
omega_thrust = outKF.omega_thrust.Data;

r = outKF.r.Data;
v = outKF.v.Data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gyroscope no noise, no bias validation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
gyro_no_noise = outKF.omega_gyro_no_noise.Data;
disp("===================")
disp("= Gyro Validation =")
disp("===================")
% Compare the gyroscope measurement model with input body rates; expected
% to be identical
err_gyro = gyro_no_noise(1,:)-w_b_ECI(1,:);
fprintf('\nInput body rates in ECI:  %e %e %e\n',w_b_ECI(1,1),w_b_ECI(1,2),w_b_ECI(1,3))
fprintf('MATLAB Three-axis gyro in ECI: %e %e %e\n',gyro_no_noise(1,1),gyro_no_noise(1,2),gyro_no_noise(1,3))
fprintf('Error vector: %e %e %e\n',err_gyro(1),err_gyro(2),err_gyro(3))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sun/Magnetometer Algorithm to calc sun_ECI before rotation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UKF_sun_t = outKF.y_UKF_sun.Time;
UKF_sun_step = nnz(UKF_sun_t==0);
% UKF_mag_t = outKF.y_UKF_mag.Time;
% UKF_mag_step = nnz(UKF_mag_t==0);
UKF_gyro_t = outKF.y_UKF_gyro.Time;
UKF_gyro_step = nnz(UKF_gyro_t==0);

time_indx_sun = ones([1 length(UKF_sun_t)/UKF_sun_step]);
% time_indx_mag = ones([1 length(UKF_mag_t)/UKF_mag_step]);
time_indx_gyro = ones([1 length(UKF_gyro_t)/UKF_gyro_step]);

count = 1;
for i = 1:UKF_sun_step:length(UKF_sun_t)
    indxs = find(UKF_sun_t(i)<=t);
    time_indx_sun(count) = indxs(1);
    count = count + 1;
end
% count = 1;
% for i = 1:UKF_mag_step:length(UKF_mag_t)
%     indxs = find(UKF_mag_t(i)<=t);
%     time_indx_mag(count) = indxs(1);
%     count = count + 1;
% end
count = 1;
for i = 1:UKF_gyro_step:length(UKF_gyro_t)
    indxs = find(UKF_gyro_t(i)<=t);
    time_indx_gyro(count) = indxs(1);
    count = count + 1;
end

% quat_mag_trunc = quat_b_ECI;
r_mag = r;
quat_sun_trunc = quat_b_ECI;
w_b_ECI_trunc = w_b_ECI(time_indx_gyro,:);



step = 100;
sun_b = zeros([floor(length(t)/step) 3]);
igrf_m_b = sun_b;
Tmb = sun_b;
T_srp = sun_b;
T_drag = sun_b;
T_a = sun_b;

grav_vec = sun_b;
g_spherical = grav_vec;
g_third_body = grav_vec;
g_tides = grav_vec;

g_schwarz = grav_vec;
g_lens = grav_vec;
g_sitter = grav_vec;
E_inertia = diag([8.0100829,8.0102594,8.0364807].*1e37); % Earth approx principal inertia tensor
w_Earth = 7.292124e-5*[0;0;1]; % rad/s earth angular velocity
EJ = (E_inertia*w_Earth)./5.97219e24; % angular momentum vector of the Earth
c = 3e8; % speed of light m/s

h_tide = 10; % max tide height in m
R_geo = 6371e3; % avg earth radius in km
mu = 398600*1000^3;
mu_sun = 1.32712440018e20; % std grav param of the sun in m^3/s^2
mu_moon = 4.9048695e12; % std grav param the moon in m^3/s^2

Cd = 2.2;               % approximate coefficent of drag
rho = 1.49e-15;        % km/m^3 US std atmosphere

for i = 1:length(sun_b)
    Ta = [0;0;0];                    % Initialzie atmo torque vector
    Tsrp = [0;0;0];                    % Initialzie atmo torque vector

    % Calculate SRP torque and measurement vectors
    JD_new = JD_start + t(i)/3600/24; % add sim time in days to start time
    jul_cent = (JD_new-2451545)/36525; % julian centuries since 2000
    lambda_S = wrapTo360(280.4606184 + 36000.77005361*jul_cent); % mean longitude of sun in deg
    M_S = wrapTo360(357.5277233 + 35999.05034*jul_cent); % mean anomaly of sun in deg
    lambda_e = wrapTo180(lambda_S + 1.914666471*sind(M_S) + 0.019994643*sind(2*M_S)); % ecliptic long of sun in deg
    
    eps = 23.439291-0.0130042*jul_cent;
    sun_ECI = [cosd(lambda_e); cosd(eps)*sind(lambda_e); sind(eps)*sind(lambda_e)];
    C_b_ECI = Functions.quat2rot(quat_b_ECI(i,:));
    sun_b(i,:) = (C_b_ECI*sun_ECI)';
    
    % Vallado SRP algorithm
    s_dir = sun_b(i,:);     % solar radiation pressure unit vector BODY frame
    for k = 1:length(dir_vector)
    dire = dir_vector{k};
        if dot(dire,s_dir) > 0
            curr_area = areas{k};
            curr_CoP_vec = CoP_vectors{k};
            for ii = 1:length(curr_area)
                Fs = -dot(dire,s_dir)*p_s*curr_area(ii).*s_dir;
                Tsrp = Tsrp + cross(curr_CoP_vec(:,ii),Fs');
            end
        end
    end
    T_srp(i,:) = Tsrp';
    
    

    % Calculate magnetic field torque and measurement vectors
    temp_t = datetime(JD_new,'ConvertFrom','juliandate');
    [ye, m, d] = ymd(temp_t);
    [h, mm, s] = hms(temp_t);
    utc_t1 = [ye m d h mm s];

    % Magnetic field effect from MATLAB
    lla = eci2lla((r_mag(i,:).*1e3),utc_t1);
    del = lla(1); % latitude
    lambda = lla(2); % longitude
    alt_g = lla(3);
    igrf_m_NED = igrfmagm(alt_g,del,lambda,decyear(ye,m,d))*1e-9;
    igrf_m_mag = norm(igrf_m_NED);
    [igrf_m_ECEF(1),igrf_m_ECEF(2),igrf_m_ECEF(3)] = ned2ecefv(igrf_m_NED(1),igrf_m_NED(2),igrf_m_NED(3),del,lambda);
    C_ECEF_ECI = dcmeci2ecef('IAU-2000/2006',utc_t1);
    igrf_m_ECI = (C_ECEF_ECI'*igrf_m_ECEF');
    igrf_m_b(i,:) = (C_b_ECI*igrf_m_ECI)';
    Tmb(i,:) = cross(m_b,igrf_m_b(i,:));

    r_sun = planetEphemeris(JD_new,'Earth','Sun').*1e3;
    r_moon = planetEphemeris(JD_new,'Earth','Moon').*1e3;
    r_ECI = r_mag(i,:);
    r_ECI_m = r_ECI.*1e3;
   
    r_ECEF = C_ECEF_ECI*(r_ECI_m)';
    r_ECEF_u = r_ECEF./norm(r_ECEF);

    % Oblation effect from MATLAB
    grav_vec(i,:) = -(mu/norm(r_ECEF)^2)*(r_ECEF/norm(r_ECEF)); % gravity vector in body frame WITHOUT spherical harmonics
    [gx, gy, gz] = gravitysphericalharmonic(r_ECEF'); % gravity vector WITH spherical harmonics
    harm_grav = [gx,gy,gz];
    g_spherical(i,:) = harm_grav-grav_vec(i,:); % Find acceleration due just to spherical changes

    % Third body effect from Yang 2021
    g_third_body(i,:) = -mu_sun*((r_ECI_m-r_sun)/norm(r_ECI_m-r_sun)^3 + r_sun/norm(r_sun)^3) -mu_moon*((r_ECI_m-r_moon)/norm(r_ECI_m-r_moon)^3 + r_moon/norm(r_moon)^3);
    
    % Tidal effect from Acedo 2016
    g_tides(i,:) = (6*mu/5/norm(r_ECEF)^4*h_tide*R_geo).*r_ECEF_u'; % maximized force from derivative of tidal potential
    
    % General Relativity formulas, from Sosnica 2021
    r_ECI_mag = norm(r_ECI_m);
     v_ECI = v(i,:).*1e3;
    v_ECI_mag = norm(v_ECI);
    g_schwarz(i,:) = mu/c^2/r_ECI_mag^3*((4*mu/r_ECI_mag - v_ECI_mag).*r_ECI_m + 4*dot(r_ECI_m,v_ECI).*v_ECI);
    g_lens(i,:) = 2*mu/c^2/r_ECI_mag^3*(3/r_ECI_mag^2*cross(r_ECI_m',v_ECI)*dot(r_ECI_m',EJ')+cross(v_ECI,EJ'));
    [R_sun,V_sun] = planetEphemeris(JD_new,'Sun','Earth');
    R_sun = R_sun.*1e3;
    V_sun = V_sun.*1e3;
    g_sitter(i,:) = cross(3*cross(-V_sun,(-mu_sun/c^2/norm(R_sun)^3*R_sun)),v_ECI);

    % Calculate drag torque and measurement vectors from Vallado
    v_b = C_b_ECI*v_ECI';
    v_dir = v_b./v_ECI_mag;     % velocity unit vector BODY frame
    Fd_coeff = -1/2*rho*v_ECI_mag^2*Cd;  % Calculate constant coeff for speed
    for p = 1:length(dir_vector)
    dire = dir_vector{p};
    if dot(dire,v_dir) >= 0
        curr_area = areas{p};
        curr_CoP_vec = CoP_vectors{p};
        for ii = 1:length(curr_area)
            Fd = dot(dire,v_dir)*Fd_coeff*curr_area(ii).*v_dir;
            Ta = Ta + cross(curr_CoP_vec(:,ii),Fd);
        end
    end
    end
    T_a(i,:) = Ta';
end

% Extract measurements from continuous timesteps corresponding to discrete steps 
y_sun = sun_b;
y_mag = igrf_m_b;
y_gyro = w_b_ECI_trunc;

% Unsure on the procedure, but Simulink EKF block calls measurements
% functions a number of times affected by measurement frequency and EKF vs
% EKF procedure. Extract set elements
y_UKF_sun_trunc = y_UKF_sun(:,1:UKF_sun_step:end);
% y_UKF_mag_trunc = y_UKF_mag(:,1:UKF_mag_step:end);
y_UKF_gyro_trunc = y_UKF_gyro(1:UKF_gyro_step:end,:);

% y_mag_no_noise = outKF.y_mag_no_noise.Data;
y_sun_no_noise = outKF.y_sun_no_noise.Data;

y_sun_no_noise = reshape(y_sun_no_noise,[3 length(y_sun_no_noise)]);
% y_mag_no_noise = reshape(y_mag_no_noise,[3 length(y_mag_no_noise)]);


err_lin_sun = zeros([length(y_sun) 3]);
for i = 1:length(y_sun)
    err_lin_sun(i,:) = y_sun(i,:)-y_sun_no_noise(:,i)';
end
% 
% err_lin_mag = zeros([length(y_mag) 3]);
% for i = 1:length(y_mag)
%     err_lin_mag(i,:) = y_mag(i,:)-y_mag_no_noise(:,i)';
% end

err_mod_gyro = zeros([length(y_gyro) 3]);
for i = 1:length(y_gyro)
    err_mod_gyro(i,:) = y_gyro(i,:)-y_UKF_gyro_trunc(i,:);
end

%% Display results
srp_torqu_norm = vecnorm(T_srp,2,2);
mag_torqu_norm = vecnorm(Tmb,2,2);
drag_torqu_norm = vecnorm(T_a,2,2);

srp_torqu_max = max(srp_torqu_norm);
mag_torqu_max = max(mag_torqu_norm);
drag_torqu_max = max(drag_torqu_norm);

grav_vec_norm = vecnorm(grav_vec,2,2);
g_sph_norm = vecnorm(g_spherical,2,2);
g_third_norm = vecnorm(g_third_body,2,2);
g_tide_norm = vecnorm(g_tides,2,2);

g_schwarz_norm = vecnorm(g_schwarz,2,2);
g_lens_norm = vecnorm(g_lens,2,2);
g_sitter_norm = vecnorm(g_sitter,2,2);

g_sph_max = max(g_sph_norm);
g_max = max(grav_vec_norm);
g_third_max = max(g_third_norm);
g_tides_max = max(g_tide_norm);

g_schwarz_max = max(g_schwarz_norm);
g_lens_max = max(g_lens_norm);
g_sitter_max = max(g_sitter_norm);

disp("=============================")
disp("= Spacecraft-Based Perturbs =")
disp("=============================")
fprintf('Max torque due to SRP: %e\n',srp_torqu_max)
fprintf('Max torque due to drag: %e\n',drag_torqu_max)
fprintf('Max torque due to mag: %e\n',mag_torqu_max)

disp("===========")
disp("= Gravity =")
disp("===========")
fprintf('\nMax grav: %e\n',g_max)
fprintf('Max grav due to oblation: %e\n',g_sph_max)
fprintf('Max grav due to third body: %e\n',g_third_max)
fprintf('Max grav due to tides: %e\n',g_tides_max)

disp("==============")
disp("= Relativity =")
disp("==============")
fprintf('\nMax Schwarzschild term: %e\n',g_schwarz_max)
fprintf('Max Lens-Thirring term: %e\n',g_lens_max)
fprintf('Max de Sitter term: %e\n',g_sitter_max)


%%
n = 1:length(y_sun);
figure
plot(n,err_lin_sun(:,1),n,err_lin_sun(:,2),n,err_lin_sun(:,3))
hold on
grid on
xlabel('n-Step')
ylabel('Sun Sensor Linearization Error')
title('Unscented Kalman Filter Sun Sensor Linearization Error')
legend('Sun X','Sun Y','Sun Z')
f = gcf;
f.Units = 'normalized';
f.Position = [0.2 0.2 0.4 0.4];
exportgraphics(gcf,'UKF_sun_line_err.png')


% n = 1:length(y_mag);
% figure
% plot(n,err_lin_mag(:,1),n,err_lin_mag(:,2),n,err_lin_mag(:,3))
% hold on
% grid on
% xlabel('n-Step')
% ylabel('Magnetometer Linearization Error')
% title('Unscented Kalman Filter Magnetometer Linearization Error')
% legend('Magnetometer X','Magnetometer Y','Magnetometer Z')
% f = gcf;
% f.Units = 'normalized';
% f.Position = [0.2 0.2 0.4 0.4];
% exportgraphics(gcf,'UKF_mag_line_err.png')

n = 1:length(y_gyro);
figure
plot(n,err_mod_gyro(:,1),n,err_mod_gyro(:,2),n,err_mod_gyro(:,3))
hold on
grid on
xlabel('n-Step')
ylabel('Gyroscope Simulation Error')
title('Unscented Kalman Filter Gyroscope Meas. Fnc Error')
legend('Gyroscope X','Gyroscope Y','Gyroscope Z')
f = gcf;
f.Units = 'normalized';
f.Position = [0.2 0.2 0.4 0.4];
exportgraphics(gcf,'UKF_gyro_line_err.png')

% figure
% plot(t,y_sun(:,1),t,y_sun(:,2),t,y_sun(:,3))
% hold on
% grid on
% plot(t,y_sun_no_noise(1,:),t,y_sun_no_noise(2,:),t,y_sun_no_noise(3,:))
% xlabel('Time (s)')
% ylabel('Sun Vector in Body Frame')
% title('Overlaid Linearized and Nonlinear Sun Vector Transformation')
% legend('Nonlin_x','Nonlin_y','Nonlin_z','Lin_x','Lin_y','Lin_z')
% f = gcf;
% f.Units = 'normalized';
% f.Position = [0.2 0.2 0.4 0.4];
% exportgraphics(gcf,'sun_combo_plot.png')

% figure
% plot(t,y_mag(:,1),t,y_mag(:,2),t,y_mag(:,3))
% hold on
% grid on
% plot(t,y_mag_no_noise(1,:),t,y_mag_no_noise(2,:),t,y_mag_no_noise(3,:))
% xlabel('Time (s)')
% ylabel('Magnetic Vector in Body Frame')
% title('Overlaid Linearized and Nonlinear Magnetic Vector Transformation')
% legend('Nonlin_x','Nonlin_y','Nonlin_z','Lin_x','Lin_y','Lin_z')
% f = gcf;
% f.Units = 'normalized';
% f.Position = [0.2 0.2 0.4 0.4];
% exportgraphics(gcf,'mag_combo_plot.png')

% figure
% plot(t1,y_gyro(:,1),t1,y_gyro(:,2),t1,y_gyro(:,3))
% hold on
% grid on
% plot(t1,y_UKF_gyro_trunc(:,1),t1,y_UKF_gyro_trunc(:,2),t1,y_UKF_gyro_trunc(:,3))
% xlabel('Time (s)')
% ylabel('Gyroscope Meas. Vector in Body Frame')
% title('Overlaid True and Measured, only Process Noise')
% legend('Nonlin_x','Nonlin_y','Nonlin_z','Lin_x','Lin_y','Lin_z')
% f = gcf;
% f.Units = 'normalized';
% f.Position = [0.2 0.2 0.4 0.4];
% exportgraphics(gcf,'gyro_combo_plot.png')
