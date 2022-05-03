function [quat] = euler2quat(eulers)
%EULER2QUAT Take a single vector (3x1, 1x3) of euler angles and convert to quaternions; expects
%to be given rotation angles as PHI, THETA, PSI, or X, Y, and Z. Calculates
%the quaternion as a 321 rotation (matrix-math wise, right-to-left).
%Extracted from MATLAB's angle2quat function.
% Jonathan Hood 2/21/2021

% Calculate half sines and cosines
cang = cos(eulers/2);
sang = sin(eulers/2);

% Calculate quaternion
% 'ZYX'
% eta = cang(3).*cang(2).*cang(1) - sang(3).*cang(2).*sang(1);
% eps =  [cang(3).*cang(2).*sang(1) + sang(3).*cang(2).*cang(1); ...
%         cang(3).*sang(2).*cang(1) + sang(3).*sang(2).*sang(1); ...
%         sang(3).*sang(2).*cang(1) - cang(3).*sang(2).*sang(1)];

% 'XYZ'
eta = cang(3).*cang(2).*cang(1) + sang(3).*sang(2).*sang(1);
eps = [cang(3).*cang(2).*sang(1) - sang(3).*sang(2).*cang(1); ...
         cang(3).*sang(2).*cang(1) + sang(3).*cang(2).*sang(1); ...
         sang(3).*cang(2).*cang(1) - cang(3).*sang(2).*sang(1)];
    
quat = [eps;eta]; % output quaternion from input euler angles
end

