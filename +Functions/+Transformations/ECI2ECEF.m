%--------------------------------------------------------------------------
%
% ECI2ECEF: Transforms Earth Centered Inertial (ECI) coordinates to Earth
%           Centered Earth Fixed (ECEF) coordinates
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function Y = ECI2ECEF(MJD_UTC, Y0)

%global eopdata const
data = load("./+Functions/+Transformations/eopdata.mat");
eopdata = data.epdata;
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = Functions.Transformations.IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = Functions.Transformations.timediff(UT1_UTC,TAI_UTC);

MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT  = MJD_UTC + TT_UTC/86400; 

% ICRS to ITRS transformation matrix and derivative
P      = Functions.Transformations.PrecMatrix(51544.5,MJD_TT);     % IAU 1976 Precession
N      = Functions.Transformations.NutMatrix(MJD_TT);                      % IAU 1980 Nutation
Theta  = Functions.Transformations.GHAMatrix(MJD_UT1);                     % Earth rotation
Pi     = Functions.Transformations.PoleMatrix(x_pole,y_pole);              % Polar motion

S = zeros(3);
S(1,2) = 1; S(2,1) = -1;                         % Derivative of Earth rotation 
Omega = 7292115.8553e-11+4.3e-15*( (MJD_UTC-51544.5)/36525 ); % [rad/s]
% Omega = const.omega_Earth-0.843994809*1e-9*LOD;  % IERS
dTheta = Omega*S*Theta; % matrix [1/s]

U      = Pi*Theta*N*P;                           % ICRS to ITRS transformation
dU     = Pi*dTheta*N*P;                          % Derivative [1/s]

% Transformation from ICRS to WGS
r = U*Y0(1:3)';
v = U*Y0(4:6)' + dU*Y0(1:3)';
Y = [r;v];

