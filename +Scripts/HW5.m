%% AERO 560 HW 5
% Jonathan Hood 3/2/2021
%% Problem 1: KF
warning('off','all')
clear
clc
close all

wn = 2;                 % natural frequency
zeta = 0.7;             % damping coefficient
% wn = 2;                 % natural frequency
% zeta = 0.6;             % damping coefficient
x0 = 0.05;              % initial position in m
xdot0 = -0.01;          % initial position in m/s
state0 = [x0; xdot0];   % set initial state
Q = 0.1;                % covariance process noise
R = 0.1;                % covariance measurement noise
variance = 0.1;
sigma_v = variance; % std. deviation of process
sigma_w = variance; % std. deviation of measurement
delta_t = 0.1;          % discrete time step in s
simtime = 10;          % simulation time in seconds
num = wn^2;
den = [1 2*wn*zeta num];

% Convert Continuous system to Discrete
[A,B,C,D] = tf2ss(num,den);
sysc = ss(A,B,C,D);
[sysd, ~] = c2d(sysc,delta_t);

A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

F_k = A;
G_k = B;
H_k = C;

L_k = G_k;
M_k = H_k;

% Initialize estimation values
x_hat0 = [0;0];
[row,col] = size(x_hat0);
P0 = eye(row);

% 1. Model Continuous time system in state space
% 2. Use c2d() command to discretize the model
% 3. Create Simulink model to simulate the MSD response for:
%   1. No noise
%   2. With noise
%   3. With noise and discrete time KF
load_system("+Sims/HW5_KF_sim"); %load simulation
set_param('HW5_KF_sim', 'StopTime', 'simtime') %set sim to stop after 1 orbital period
outKF = sim('HW5_KF_sim');      % Run Simulink model

% Plot MSD for all three systems
noiseNoFilter = outKF.noiseNoFilter.Data;
noNoise = outKF.noNoise.Data;
x_hat_k = outKF.x_hat_k.Data;
y_hat_k = outKF.y_hat_k.Data;
t = 0:delta_t:simtime;

figure
plot(t,noiseNoFilter)
hold on
grid on
xlabel('Time (s)')
ylabel('Signal')
title('Signal Filter Comparison')
plot(t,noNoise)
plot(t,y_hat_k(:,1))
legend('Noise w/o KF','w/o Noise, w/o KF','Noise w/ KF')

%% Problem 2: Norm-Constrained EKF

clear
close all
clc

%%%%%%%%%%%%%%
% S/C Givens %
%%%%%%%%%%%%%%

% Uncontrolled s/c
u = 0; % for all k, no control signal

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

%%%%%%%%%%%%%%%%%%%%%%%
% Calculate rotations %
%%%%%%%%%%%%%%%%%%%%%%%

% Initial b_ECI
w_b_ECI = [0.05; -0.05; 0.05]; % initial angular velocity rad/s
quat_b_ECI = [sin(1/2)/sqrt(3)*ones([3 1]); cos(1/2)]; % initial quaternion
C_b_ECI = Functions.quat2rot(quat_b_ECI); % body ECI rotation matrix
euler_b_ECI = Functions.rot2euler(C_b_ECI); % body ECI euler angles

% %%%%%%%%%%%%% CALCULATE INITIAL ROTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
z = -r/norm(r);                                % z basis vector LVLH
y = -cross(r,v)/norm(cross(r,v));              % y basis vector LVLH
x = cross(y,z);                                % x basis vector LVLH

C_LVLH_ECI = [x;y;z];                             % get ECI to LVLH rotation matrix 
quat_LVLH_ECI = Functions.rot2quat(C_LVLH_ECI);   % get ECI to LVLH quaternion
euler_LVLH_ECI = Functions.rot2euler(C_LVLH_ECI); % get ECI to LVLH euler angles

% Initial b_LVLH
C_b_LVLH = C_b_ECI*C_LVLH_ECI';                     % get b to LVLH rotation matrix 
quat_b_LVLH = Functions.rot2quat(C_b_LVLH);         % get b to LVLH quaternion
euler_b_LVLH = Functions.rot2euler(C_b_LVLH);       % get b to LVLH euler angles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Given Sensor/Kalman Conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Norm-Constrained EKF initial conditions
w_hat0 = [0;0;0]; % initial angular velocity
quat_hat0 = [0;0;0;1]; % initial quaternion
x_hat0 = [w_hat0; quat_hat0]; % initial state
P0 = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.1]); % initial covariance

% Process Noise
sigma_q = 0.001; % std deviation process noise
var_q = sigma_q^2; % process noise variance
Q_k = sigma_q^2*eye(); % process noise covariance

% Measurement noise
sigma_m = 0.01; % std deviation magnetometer measurement noise
var_m = sigma_m^2; % magnetometer measurement noise variance
sigma_s = 0.005; % std deviation sun sensor measurement noise
var_s = sigma_s^2; % sun sensor measurement noise variance
R_k = diag([sigma_m^2,sigma_m^2,sigma_m^2,sigma_s^2,sigma_s^2,sigma_s^2]); % initial measurement noise covariance

M_k = eye(size(R_k)); % identity matrix size of R_k

delta_t = 1; % time step between measurements & prediction time step

L_k = [delta_t*inv(J); zeros(3); zeros([1 3])];

% Sun and Magnetometer Measurement vectors
sun_ECI = [1;0;0]; % sun measurement vector in ECI
sun_hat_ECI = sun_ECI./norm(sun_ECI);
mag_ECI = [0.3815; -0.0969; 0.9193]; % magnetic field measurement vector in ECI
mag_hat_ECI = mag_ECI./norm(mag_ECI);

load_system("+Sims/HW5_NCEKF_sim"); %load simulation
set_param('HW5_NCEKF_sim', 'StopTime', 'simtime') %set sim to stop after 1 orbital period
outKF = sim('HW5_NCEKF_sim');      % Run Simulink model

% Plot MSD for all three systems
w_b_ECI = outKF.w_b_ECI.Data;
x_hat_k = outKF.x_hat_k.Data;
w_b_ECI_filt = x_hat_k(:,1:3);
quat_b_ECI_filt = x_hat_k(:,4:7);
quat_b_ECI = outKF.quats_b_ECI.Data;
t = outKF.tout;

% Plot angular velocity
figure
plot(t,w_b_ECI(:,1),t,w_b_ECI(:,2),t,w_b_ECI(:,3))
hold on
grid on
xlabel('Time (s)')
ylabel('\omega (rad/s)')
title('Angular Velocity: Non-filtered vs. Filtered')
plot(t,w_b_ECI_filt(:,1),'--',t,w_b_ECI_filt(:,2),'--',...
    t,w_b_ECI_filt(:,3),'--')
legend('\omega_1','\omega_2','\omega_3','\omega_{F1}','\omega_{F2}',...
    '\omega_{F3}')

% Plot quaternion components
figure
plot(t,quat_b_ECI(:,1),t,quat_b_ECI(:,2),t,quat_b_ECI(:,3),...
    t,quat_b_ECI(:,4))
hold on
grid on
xlabel('Time (s)')
ylabel('Magnitude')
title('Quaternion: Non-filtered vs. Filtered')
plot(t,quat_b_ECI_filt(:,1),'--',t,quat_b_ECI_filt(:,2),'--',...
    t,quat_b_ECI_filt(:,3),'--',t,quat_b_ECI_filt(:,4),'--')
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta','\epsilon_{F1}',...
    '\epsilon_{F2}','\epsilon_{F3}','\eta_F')



