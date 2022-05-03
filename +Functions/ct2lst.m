function lst = ct2lst(lon, JD)
% CT2LST  Find the local sidereal time from the local Civilian Time in
% degrees and Julian Date
%
% lst = ct2lst(lon, JD)
% Returns the local sidereal time in degrees, lst, given the local civilian
% time in degrees, lon, and the Julian Date, JD.  The local Civilian Time
% is the same as the East Longitude of the current location of interest on
% earth.
%
% From Curtis
% Author: Eric A. Mehiel
% Date: December 19th, 2006

J0 = floor(JD + 0.5) - 0.5;
UT = (JD - J0)*24.0;
T0 = (J0 - 2451545.0)/36525.0;
c = [100.4606184, 36000.77004, 0.000387933, -2.583e-8, 360.98564724];
g0 = c(1) + c(2)*T0 + c(3)*T0^2 + c(4)*T0^3;

g0 = g0 - 360.0*floor(g0/360.0);

gst = g0 + c(5)*UT/24.0;
lst = gst + lon;

lst = lst - 360.0*floor(lst/360.0);