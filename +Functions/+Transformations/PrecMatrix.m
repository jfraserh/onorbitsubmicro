%--------------------------------------------------------------------------
%
% PrecMatrix: Precession transformation of equatorial coordinates
%
% Inputs:
%   Mjd_1     Epoch given (Modified Julian Date TT)
%   MjD_2     Epoch to precess to (Modified Julian Date TT)
% 
% Output:
%   PrecMat   Precession transformation matrix
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function PrecMat = PrecMatrix(Mjd_1, Mjd_2)


% Mathematical constants
const.pi2       = 2*pi;                % 2pi
const.Rad       = pi/180;              % Radians per degree
const.Deg       = 180/pi;              % Degrees per radian
const.Arcs      = 3600*180/pi;         % Arcseconds per radian

% General
const.MJD_J2000 = 51544.5;             % Modified Julian Date of J2000
const.T_B1950   = -0.500002108;        % Epoch B1950
const.c_light   = 299792458.000000;    % Speed of light  [m/s]; DE430
const.AU        = 149597870700.000000; % Astronomical unit [m]; DE430

% Physical parameters of the Earth, Sun and Moon

% Equatorial radius and flattening
const.R_Earth   = 6378.137e3;          % Earth's radius [m]; WGS-84
const.f_Earth   = 1/298.257223563;     % Flattening; WGS-84   
const.R_Sun     = 696000e3;            % Sun's radius [m]; DE430
const.R_Moon    = 1738e3;              % Moon's radius [m]; DE430

% Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
const.omega_Earth = 15.04106717866910/3600*pi/180; % [rad/s]; WGS-84

% Gravitational coefficients
const.GM_Earth   = 398600.4418e9;               	 % [m^3/s^2]; WGS-84
const.GM_Sun     = 132712440041.939400e9;      	  	 % [m^3/s^2]; DE430
const.GM_Moon    = 398600.4418e9/81.30056907419062; % [m^3/s^2]; DE430
const.GM_Mercury = 22031.780000e9;             	  	 % [m^3/s^2]; DE430
const.GM_Venus   = 324858.592000e9;            	  	 % [m^3/s^2]; DE430
const.GM_Mars    = 42828.375214e9;             	  	 % [m^3/s^2]; DE430
const.GM_Jupiter = 126712764.800000e9;         	  	 % [m^3/s^2]; DE430
const.GM_Saturn  = 37940585.200000e9;          	  	 % [m^3/s^2]; DE430
const.GM_Uranus  = 5794548.600000e9;           	  	 % [m^3/s^2]; DE430
const.GM_Neptune = 6836527.100580e9;           	  	 % [m^3/s^2]; DE430
const.GM_Pluto   = 977.0000009e9;        	  	 	 % [m^3/s^2]; DE430

% Solar radiation pressure at 1 AU 
const.P_Sol = 1367/299792458.000000; % [N/m^2] (1367 W/m^2); IERS 96



T  = (Mjd_1-const.MJD_J2000)/36525;
dT = (Mjd_2-Mjd_1)/36525;

% Precession angles
zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+ ...
           ((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/const.Arcs;
z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/const.Arcs;
theta =  ( (2004.3109-(0.85330+0.000217*T)*T)- ...
           ((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/const.Arcs;

% Precession matrix
PrecMat = Functions.Transformations.R_z(-z) * Functions.Transformations.R_y(theta) * Functions.Transformations.R_z(-zeta);

