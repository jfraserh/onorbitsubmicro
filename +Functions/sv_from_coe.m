% 
function [r, v] = sv_from_coe(coe,mu)
% 
%{
This function computes the state vector (r,v) from the
classical orbital elements (coe).
mu - gravitational parameter (km^3/s^2)
coe - orbital elements [h e RA incl w TA]
where
h = angular momentum (km^2/s)
e = eccentricity
RA = right ascension of the ascending node (deg)
incl = inclination of the orbit (deg)
w = argument of perigee (deg)
TA = true anomaly (deg)
R3_w - Rotation matrix about the z-axis through the angle w
R1_i - Rotation matrix about the x-axis through the angle i
R3_W - Rotation matrix about the z-axis through the angle RA
Q_pX - Matrix of the transformation from perifocal to geocentric
equatorial frame
rp - position vector in the perifocal frame (km)
vp - velocity vector in the perifocal frame (km/s)
r - position vector in the geocentric equatorial frame (km)
v - velocity vector in the geocentric equatorial frame (km/s)
MATLAB scripts e59
User M-functions required: none
%}
% ––––––––––––––––––––––––––––––––––––––––––––––
h = coe(1);
e = coe(2);
RA = coe(3);
incl = coe(4);
w = coe(5);
TA = coe(6);
%...Equations 4.45 and 4.46 (rp and vp are column vectors):
rp = (h^2/mu) * (1/(1 + e*cosd(TA))) * (cosd(TA)*[1;0;0] + sind(TA)*[0;1;0]);
vp = (mu/h) * (-sind(TA)*[1;0;0] + (e + cosd(TA))*[0;1;0]);
%...Equation 4.34:
R3_W = [ cosd(RA) sind(RA) 0
-sind(RA) cosd(RA) 0
0 0 1];
%...Equation 4.32:
R1_i = [1 0 0
0 cosd(incl) sind(incl)
0 -sind(incl) cosd(incl)];
%...Equation 4.34:
R3_w = [ cosd(w) sind(w) 0
-sind(w) cosd(w) 0
0 0 1];
%...Equation 4.49:
Q_pX = (R3_w*R1_i*R3_W)';
%...Equations 4.51 (r and v are column vectors):
r = Q_pX*rp;
v = Q_pX*vp;
%...Convert r and v into row vectors:
r = r';
v = v';
end
% 