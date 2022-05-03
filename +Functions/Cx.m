function [C_x] = Cx(fi)
%CX Summary of this function goes here
%   Detailed explanation goes here
C_x =[1  0       0;
      0  cos(fi) sin(fi);
      0 -sin(fi) cos(fi)];
end

