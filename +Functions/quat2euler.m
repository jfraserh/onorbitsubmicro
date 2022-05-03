function [euler] = quat2euler(quaternion)
%Turn input quaternion into euler angles
% Reference:
% https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Quaternion_to_Euler_Angles_Conversion
q0 = quaternion(4);
q1 = quaternion(1);
q2 = quaternion(2);
q3 = quaternion(3);

fi = atan2(2*(q0*q1+q2*q3),1-2*(q1^2+q2^2));
theta = asin(2*(q0*q2-q3*q1));
psi = atan2(2*(q0*q3+q1*q2),1-2*(q2^2+q3^2));

euler = [fi;theta;psi];

end

