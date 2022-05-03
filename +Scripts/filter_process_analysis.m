%% Filter Process and Measurement Noise Tradespace
% Load all process and measurement noise runs and save minimum at UKF and
% EKF. Display global minimum related initial covariance, process noise,
% and related sun sensor and gyroscope measurement noise. Run
% Filter_Adjust_Microthruster_Script.m to generate the appropriate data.

clear
clc 
close all
EKF_arr = [];
UKF_arr = [];
str_arr = ["08","07","06","05","0001","001","01","1"];

q_arr = ["08","07","06","05","0001","001","01","1"];
countttt = 0;
for iii = 1:length(str_arr)
    for jjj = 1:length(q_arr)
        load(strjoin(['Filter_P0_' str_arr(iii) '_Meas_all_' q_arr(jjj) '_outKF_Full_Sim.mat'],''));
        
        [minv,indx] =  min(abs(thrust_mean_EKF.*1e6-0.1),[],'all'); % linear index
        EKF_arr = [EKF_arr; minv indx];
        [minv,indx] =  min(abs(thrust_mean_UKF.*1e6-0.1),[],'all'); % linear index
        UKF_arr = [UKF_arr; minv indx];
    end
end

%% Find mins and respective Q and P0 values
[EKF_min,E_indxs] = min(EKF_arr(:,1));
[UKF_min,U_indxs] = min(UKF_arr(:,1));

Q_val = str_arr(mod(E_indxs,length(q_arr)));
P0_val = str_arr(ceil(E_indxs/length(str_arr)));

disp("==========")
disp("EKF Result")
disp("==========")

fprintf("\nProcess Q: %s",Q_val)
fprintf("\nInitial Cov.: %s",P0_val)
EKF_lin_indx = EKF_arr(E_indxs,2);
fprintf("\nThrust value-0.1, index: %f, %d",EKF_min,EKF_lin_indx)

% Measurement noise matrices is a 5x5; linear index, convert to row and
% column
col = ceil(EKF_lin_indx/5); % Measurement matrix row
row = mod(EKF_lin_indx,5); % Measurement matrix col
gyro_opt_noise = meas_rang_gyro(row);
sun_opt_noise = meas_rang_sun(col);

fprintf("\nGyro measurement noise: %f",gyro_opt_noise)
fprintf("\nSun measurement noise: %f\n\n",sun_opt_noise)

disp("==========")
disp("UKF Result")
disp("==========")

fprintf("\nProcess Q: %s",Q_val)
fprintf("\nInitial Cov.: %s",P0_val)
UKF_lin_indx = UKF_arr(E_indxs,2);
fprintf("\nThrust value-0.1, index: %f, %d",UKF_min,UKF_lin_indx)

% Measurement noise matrices is a 5x5; linear index, convert to row and
% column
col = ceil(UKF_lin_indx/5); % Measurement matrix row
row = mod(UKF_lin_indx,5); % Measurement matrix col
gyro_opt_noise = meas_rang_gyro(row);
sun_opt_noise = meas_rang_sun(col);

fprintf("\nGyro measurement noise: %f",gyro_opt_noise)
fprintf("\nSun measurement noise: %f\n\n",sun_opt_noise)
