%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%  
% Last modified:   2018/02/01   M. Mahooti
% 
%--------------------------------------------------------------------------
function [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eop,Mjd_UTC,interp)

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


if (nargin == 2)
   interp = 'n';
end

if (interp =='l')
    % linear interpolation
    mjd = (floor(Mjd_UTC));
    i = find(mjd==eop(:,1),1,'first');
    preeop = eop(i,:);
    nexteop = eop(i+1,:);
    mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
    fixf = mfme/1440;
    % Setting of IERS Earth rotation parameters
    % (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    x_pole  = preeop(2)+(nexteop(2)-preeop(2))*fixf;
    y_pole  = preeop(3)+(nexteop(3)-preeop(3))*fixf;
	UT1_UTC = preeop(4)+(nexteop(4)-preeop(4))*fixf;
    LOD     = preeop(5)+(nexteop(5)-preeop(5))*fixf;
    dpsi    = preeop(6)+(nexteop(6)-preeop(6))*fixf;
    deps    = preeop(7)+(nexteop(7)-preeop(7))*fixf;
    dx_pole = preeop(8)+(nexteop(8)-preeop(8))*fixf;
    dy_pole = preeop(9)+(nexteop(9)-preeop(9))*fixf;
    TAI_UTC = preeop(10);
	
    x_pole  = x_pole/const.Arcs;  % Pole coordinate [rad]
    y_pole  = y_pole/const.Arcs;  % Pole coordinate [rad]
    dpsi    = dpsi/const.Arcs;
    deps    = deps/const.Arcs;
    dx_pole = dx_pole/const.Arcs; % Pole coordinate [rad]
    dy_pole = dy_pole/const.Arcs; % Pole coordinate [rad]
elseif (interp =='n')    
    mjd = (floor(Mjd_UTC));
    i = find(mjd==eop.MJD,1,'first');
    eop = eop(:,i);
    % Setting of IERS Earth rotation parameters
    % (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    x_pole  = eop(5)/const.Arcs;  % Pole coordinate [rad]
    y_pole  = eop(6)/const.Arcs;  % Pole coordinate [rad]
	UT1_UTC = eop(7);             % UT1-UTC time difference [s]
    LOD     = eop(8);             % Length of day [s]
    dpsi    = eop(9)/const.Arcs;
    deps    = eop(10)/const.Arcs;
    dx_pole = eop(11)/const.Arcs; % Pole coordinate [rad]
    dy_pole = eop(12)/const.Arcs; % Pole coordinate [rad]
	TAI_UTC = eop(13);            % TAI-UTC time difference [s]
end

