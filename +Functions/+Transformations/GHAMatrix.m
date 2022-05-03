%--------------------------------------------------------------------------
%
% GHAMatrix: Transformation from true equator and equinox to Earth equator
%            and Greenwich meridian system 
%
% Input:
%   Mjd_UT1   Modified Julian Date UT1
% 
% Output:
%   GHAmat    Greenwich Hour Angle matrix
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function GHAmat = GHAMatrix(Mjd_UT1)

GHAmat = Functions.Transformations.R_z(Functions.Transformations.gast(Mjd_UT1));

