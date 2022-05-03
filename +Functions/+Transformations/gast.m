%--------------------------------------------------------------------------
%
% gast: Greenwich Apparent Sidereal Time
%
% Input:
%   Mjd_UT1   Modified Julian Date UT1
%
% Output:
%   gstime    GAST in [rad]
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function gstime = gast(Mjd_UT1)

gstime = mod(Functions.Transformations.gmst(Mjd_UT1) + Functions.Transformations.EqnEquinox(Mjd_UT1), 2*pi);

