function [cross] = vec2cross(vec)
%VEC2CROSS Summary of this function goes here
%   Detailed explanation goes here
cross = [0    -vec(3)   vec(2);
        vec(3) 0       -vec(1);
       -vec(2) vec(1)   0];
end

