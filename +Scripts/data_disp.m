%% Data Display Script
% Displays datasets from filter process runs.

clear
clc 
close all

EKF_arr = [];
UKF_arr = [];
str_arr = ["08","07","06","05","0001","001","01","1"];
num_arr = [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];
q_arr = ["08","07","06","05","0001","001","01","1"];

for iii = 1:length(str_arr)
    for jjj = 1:length(q_arr)
        load(strjoin(['Filter_P0_' str_arr(iii) '_Meas_all_' q_arr(jjj) '_outKF_Full_Sim.mat'],''));
        fprintf('\nP0 = %.2e, Q = %.2e \n',num_arr(iii),num_arr(jjj))
        disp("===========================")
        disp("EKF")
        disp(thrust_mean_EKF.*1e6)
        disp("UKF")
        disp(thrust_mean_UKF.*1e6)
    end
end