function [quat] = rot2quat(rot)
%rot2quat Summary of this function goes here
%   Detailed explanation goes here
eta = sqrt(trace(rot) + 1)/2;
tol = 2e-6;
if eta~=0
    eps1 = (rot(2,3)-rot(3,2))/(4*eta);
    eps2 = (rot(3,1)-rot(1,3))/(4*eta);
    eps3 = (rot(1,2)-rot(2,1))/(4*eta);
else
    eps1 = ((rot(1,1) + 1)/2)^(1/2);
    eps2 = ((rot(2,2) + 1)/2)^(1/2);
    eps3 = ((rot(3,3) + 1)/2)^(1/2);
    if eps1 < tol
        eps2 = sign(rot(1,2))*eps2;
        eps3 = sign(rot(3,1))*eps3;
    elseif eps2 < tol
        eps1 = sign(rot(1,2))*eps1;
        eps3 = sign(rot(3,2))*eps3;
    elseif eps3 < tol
        eps2 = sign(rot(2,3))*eps2;
        eps1 = sign(rot(3,1))*eps1;
    end
end
quat = [eps1; eps2; eps3; eta];
% take negatives into account?
%{
m00 = rot(1,1);
m01 = rot(1,2);
m02 = rot(1,3);

m10 = rot(2,1);
m11 = rot(2,2);
m12 = rot(2,3);

m20 = rot(3,1);
m21 = rot(3,2);
m22 = rot(3,3);

 tr = m00 + m11 + m22;

if (tr > 0)
   S = sqrt(tr+1.0) * 2; 
  qw = 0.25 * S;
  qx = (m21 - m12) / S;
  qy = (m02 - m20) / S; 
  qz = (m10 - m01) / S; 
elseif ((m00 > m11) && (m00 > m22)) 
   S = sqrt(1.0 + m00 - m11 - m22) * 2; 
  qw = (m21 - m12) / S;
  qx = 0.25 * S;
  qy = (m01 + m10) / S; 
  qz = (m02 + m20) / S; 
elseif (m11 > m22)
   S = sqrt(1.0 + m11 - m00 - m22) * 2; 
  qw = (m02 - m20) / S;
  qx = (m01 + m10) / S; 
  qy = 0.25 * S;
  qz = (m12 + m21) / S; 
 else 
   S = sqrt(1.0 + m22 - m00 - m11) * 2; 
  qw = (m10 - m01) / S;
  qx = (m02 + m20) / S;
  qy = (m12 + m21) / S;
  qz = 0.25 * S;
end
qtot = [qx; qy; qz; qw];
%}
end

