function [C] = quat2rot(q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
qr = q(4);
% a = qr;
qi = q(1);
% b = qi;
qj = q(2);
% c = qj;
qk = q(3);
% d = qk;
s = norm(q)^-2;

C = [1-2*s*(qj^2+qk^2) 2*s*(qi*qj-qk*qr) 2*s*(qi*qk+qj*qr);
    2*s*(qi*qj+qk*qr)  1-2*s*(qi^2+qk^2) 2*s*(qj*qk-qi*qr);
    2*s*(qi*qk-qj*qr)  2*s*(qj*qk+qi*qr) 1-2*s*(qi^2+qj^2)]';

% C1 = [2*a^2-1+2*b^2    2*b*c+2*a*d      2*b*d-2*a*c;
%       2*b*c-2*a*d      2*a^2-1+2*c^2    2*c*d+2*a*b;
%       2*b*d+2*a*c      2*c*d-2*a*b      2*a^2-1+2*d^2];

end

