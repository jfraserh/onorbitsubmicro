function [C_y] = Cy(theta)
%Cy Summary of this function goes here
%   Detailed eyplanation goes here
C_y =[cos(theta)  0         -sin(theta);
      0           1         0;
      sin(theta)  0         cos(theta)];
end