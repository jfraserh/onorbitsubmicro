function [C_z] = Cz(psi)
%Cz Summarz of this function goes here
%   Detailed ezplanation goes here
C_z =[cos(psi)  sin(psi)  0;
     -sin(psi)  cos(psi)  0;
      0           0       1];
end