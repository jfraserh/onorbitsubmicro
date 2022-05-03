function [euler] = rot2euler(rot)
%rot2euler Summary of this function goes here
%   Detailed explanation goes here
fi = atan2(rot(2,3),rot(3,3));
theta = asin(-rot(1,3));
psi = atan2(rot(1,2),rot(1,1));
euler = [fi; theta; psi];
end

