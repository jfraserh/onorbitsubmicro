%--------------------------------------------------------------------------
%
% NutMatrix: Transformation from mean to true equator and equinox
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   NutMat    Nutation matrix
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function NutMat = NutMatrix(Mjd_TT)

% Mean obliquity of the ecliptic
ep = Functions.Transformations.MeanObliquity(Mjd_TT);

% Nutation in longitude and obliquity
[dpsi, deps] = Functions.Transformations.NutAngles(Mjd_TT);

% Transformation from mean to true equator and equinox
NutMat = Functions.Transformations.R_x(-ep-deps)*Functions.Transformations.R_z(-dpsi)*Functions.Transformations.R_x(ep);

