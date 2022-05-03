%% Uncertainty Analysis
% Load uncertainty data from a series of simulation runs at different run
% times and fit a square-root curve to uncertainty. Calculate how many data
% points, and therefore how long a simtime, would be required to breach
% minimum uncertainty from 0.09 to 0.012 uN.
% Requires data from "Optimal_Microthruster_Script" in simulation test
% mode.
%

clear
clc
close all

load Optimal_Long_Runs_outKF_Full_Sim.mat

%% EKF

EKF_uncert_fit = fit(simtime_range',(thrust_uncert_EKF+2)','power1');
% Power-1 refers to a*x^b; expecting an approximately square root fit

coeffs_EKF = coeffvalues(EKF_uncert_fit);
ci_EKF = confint(EKF_uncert_fit,0.68);
a_EKF = coeffs_EKF(1);
b_EKF = coeffs_EKF(2);
% c_EKF = coeffs_EKF(3);
EKF_a_unc = abs(a_EKF-ci_EKF(1,1)); % exponent deviation
EKF_b_unc = abs(b_EKF-ci_EKF(1,2)); % exponent deviation
% EKF_c_unc = abs(c_EKF-ci_EKF(1,3)); % exponent deviation

% fprintf('\nEKF Uncertainty Equation:\n\t %fx^(%f) + %f \n',...
%     a_EKF,b_EKF,c_EKF)

fprintf('\nEKF Uncertainty Equation:\n\t %fx^(%f)\n',...
    a_EKF,b_EKF)

min_unc = 0.09;

% EKF_req_time = ((min_unc-c_EKF)/a_EKF)^(1/b_EKF);
EKF_req_time = (min_unc/a_EKF)^(1/b_EKF);

fprintf('\nRequired EKF simulation time: %f hrs, %f days\n',EKF_req_time/3600,EKF_req_time/3600/24)
fprintf('Corresponding to a 50 Hz, gyroscope, or %f data points required\n',50*EKF_req_time)

%% UKF

UKF_uncert_fit = fit(simtime_range',thrust_uncert_UKF','power1');
% Power-2 refers to a*x^b+c; expecting an approximately square root fit

coeffs_UKF = coeffvalues(UKF_uncert_fit);
ci_UKF = confint(UKF_uncert_fit,0.68);
a_UKF = coeffs_UKF(1);
b_UKF = coeffs_UKF(2);
% c_UKF = coeffs_UKF(3);
UKF_a_unc = abs(a_UKF-ci_UKF(1,1)); % exponent deviation
UKF_b_unc = abs(b_UKF-ci_UKF(1,2)); % exponent deviation
% UKF_c_unc = abs(c_UKF-ci_UKF(1,3)); % exponent deviation

% fprintf('\nUKF Uncertainty Equation:\n\t %fx^(%f) + %f \n',...
%     a_UKF,b_UKF,c_UKF)

fprintf('\nUKF Uncertainty Equation:\n\t %fx^(%f)\n',...
    a_UKF,b_UKF)


% UKF_req_time = ((min_unc-c_UKF)/a_UKF)^(1/b_UKF);
UKF_req_time = (min_unc/a_UKF)^(1/b_UKF);


fprintf('\nRequired UKF simulation time: %f s, %f min, %f hrs\n',UKF_req_time,UKF_req_time/60,UKF_req_time/3600)
fprintf('Corresponding to a 50 Hz, gyroscope, or %f data points required\n',50*UKF_req_time)
x = 347005:10000:19008000;
% ekf_time = a_EKF.*x.^b_EKF + c_EKF;
% ukf_time = a_UKF.*x.^b_UKF + c_UKF;

ekf_time = a_EKF.*x.^b_EKF;
ukf_time = a_UKF.*x.^b_UKF;
figure
loglog(x,ekf_time)
grid on
hold on
loglog(x,ukf_time)
xlabel('Simulation Time (s)')
ylabel('Thrust Uncertainty uN')
title('Thrust Uncertainty over Time')
legend('EKF','UKF')
ylim([0 0.1])
xlim([-inf inf])
exportgraphics(gcf,'thrust_uncert_graph.png')