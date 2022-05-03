function [Y_b] = Yb_func(quat,vec)
%CX Summary of this function goes here
%   Detailed explanation goes here
eps = quat(1:3);
eta = quat(4);
Y_b = [eta*eye(3)-Functions.vec2cross(eps) -eps]*[Functions.vec2cross(vec) vec;
                                                      -vec'                    0];

end

