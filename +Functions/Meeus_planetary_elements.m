function [planet_coes] = Meeus_planetary_elements(planet_id,T)
% Planetary Ephemerides from Meeus (1991:202-204) and J2000.0
% Output:
% planet_coes
% a = semimajor axis (km)
% ecc = eccentricity
% inc = inclination (degrees)
% raan = right ascension of the ascending node (degrees)
% w_hat = longitude of perihelion (degrees)
% L = mean longitude (degrees)
%
% Inputs:
% planet_id - planet identifier:
% 1 = Mercury
% 2 = Venus
% 3 = Earth
% 4 = Mars
% 5 = Jupiter
% 6 = Saturn
% 7 = Uranus
% 8 = Neptune
%
% T = Julian date in days

% Turn Julian days into Julian centuries
T = (T - 2451545.0)/36525;

% Set grav param of the sun
uS = 132712440018;

if planet_id == 1 % Mercury
    planet_coes.u = 22032; % km^3/s std. grav param
    planet_coes.a = 0.387098310; % AU but in km later
    planet_coes.ecc = 0.20563175 + 0.000020406*T - 0.0000000284*T^2 - 0.00000000017*T^3;
    planet_coes.inc = 7.004986 - 0.0059516*T + 0.00000081*T^2 + 0.000000041*T^3; %degs
    planet_coes.raan = 48.330893 - 0.1254229*T-0.00008833*T^2 - 0.000000196*T^3; %degs
    planet_coes.w_hat = 77.456119 +0.1588643*T -0.00001343*T^2+0.000000039*T^3; %degs
    planet_coes.L = 252.250906+149472.6746358*T-0.00000535*T^2+0.000000002*T^3; %degs
elseif planet_id == 2 % Venus
    planet_coes.u = 324859; % km^3/s std. grav param
    planet_coes.a = 0.723329820; % AU
    planet_coes.ecc = 0.00677188 - 0.000047766*T + 0.000000097*T^2 + 0.00000000044*T^3;
    planet_coes.inc = 3.394662 - 0.0008568*T - 0.00003244*T^2 + 0.000000010*T^3; %degs
    planet_coes.raan = 76.679920 - 0.2780080*T-0.00014256*T^2 - 0.000000198*T^3; %degs
    planet_coes.w_hat = 131.563707 +0.0048646*T -0.00138232*T^2-0.000005332*T^3; %degs
    planet_coes.L = 181.979801+58517.8156760*T+0.00000165*T^2-0.000000002*T^3; %degs
elseif planet_id == 3  % Earth
    planet_coes.u = 398600; % km^3/s std. grav param
    planet_coes.a = 1.000001018; % AU
    planet_coes.ecc = 0.01670862 - 0.000042037*T - 0.0000001236*T^2 + 0.00000000004*T^3;
    planet_coes.inc = 0.0000000 + 0.0130546*T - 0.00000931*T^2 - 0.000000034*T^3; %degs
    planet_coes.raan = 0.0; %degs
    planet_coes.w_hat = 102.937348 + 0.3225557*T + 0.00015026*T^2 + 0.000000478*T^3; %degs
    planet_coes.L = 100.466449 + 35999.372851*T - 0.00000568*T^2 + 0.000000000*T^3; %degs
elseif planet_id == 4 % Mars
    planet_coes.u = 4.2828e+04; % km^3/s std. grav param
    planet_coes.a = 1.523679342; % AU
    planet_coes.ecc = 0.09340062 + 0.000090483*T - 0.00000000806*T^2 - 0.00000000035*T^3;
    planet_coes.inc = 1.849726 - 0.0081479*T - 0.00002255*T^2 - 0.000000027*T^3; %degs
    planet_coes.raan = 49.558093 - 0.2949846*T-0.00063993*T^2 - 0.000002143*T^3; %degs
    planet_coes.w_hat = 336.060234 +0.4438898*T -0.00017321*T^2+0.000000300*T^3; %degs
    planet_coes.L = 355.433275+19140.2993313*T+0.00000261*T^2-0.000000003*T^3; %degs
elseif planet_id == 5 % Jupiter
    planet_coes.u =  126686500; % km^3/s std. grav param
    planet_coes.a = 5.202603191 + 0.0000001913*T; % AU
    planet_coes.ecc = 0.04849485+0.000163244*T - 0.0000004719*T^2 + 0.00000000197*T^3;
    planet_coes.inc = 1.303270 - 0.0019872*T + 0.00003318*T^2 + 0.000000092*T^3; %degs
    planet_coes.raan = 100.464441 + 0.1766828*T+0.00090387*T^2 - 0.000007032*T^3; %degs
    planet_coes.w_hat = 14.331309 +0.2155525*T +0.00072252*T^2-0.000004590*T^3; %degs
    planet_coes.L = 34.351484+3034.9056746*T-0.00008501*T^2+0.000000004*T^3; %degs
elseif planet_id == 6 % Saturn
    planet_coes.u =  37931190; % km^3/s std. grav param
    planet_coes.a = 9.5549009596 - 0.0000021389*T; % AU
    planet_coes.ecc = 0.05550862 - 0.000346818*T -0.0000006456*T^2 + 0.00000000338*T^3;
    planet_coes.inc = 2.488878 + 0.0025515*T - 0.00004903*T^2 + 0.000000018*T^3; %degs
    planet_coes.raan = 113.665524 - 0.2566649*T-0.00018345*T^2 + 0.000000357*T^3; %degs
    planet_coes.w_hat = 93.056787 +0.5665496*T +0.00052809*T^2-0.000004882*T^3; %degs
    planet_coes.L = 50.077471+1222.1137943*T+0.00021004*T^2-0.000000019*T^3; %degs
elseif planet_id == 7 % Uranus
    planet_coes.u =  5793940; % km^3/s std. grav param
    planet_coes.a = 19.218446062-0.0000000372*T+0.00000000098*T^2; % AU
    planet_coes.ecc = 0.04629590 - 0.000027337*T + 0.0000000790*T^2 + 0.00000000025*T^3;
    planet_coes.inc = 0.773196 - 0.0016869*T + 0.00000349*T^2 + 0.00000000016*T^3; %degs
    planet_coes.raan = 74.005947 + 0.0741461*T+0.00040540*T^2 +0.000000104*T^3; %degs
    planet_coes.w_hat = 173.005159 +0.0893206*T -0.00009470*T^2+0.000000413*T^3; %degs
    planet_coes.L = 314.055005+428.4669983*T-0.00000486*T^2-0.000000006*T^3; %degs
elseif planet_id == 8 % Neptune
    planet_coes.u =  6836520; % km^3/s std. grav param
    planet_coes.a = 30.110386869-0.0000001663*T+0.00000000069*T^2; % AU
    planet_coes.ecc = 0.00898809 + 0.000006408*T -0.0000000008*T^2;
    planet_coes.inc = 1.769952 +0.0002557*T +0.00000023*T^2 -0.0000000000*T^3; %degs
    planet_coes.raan = 131.784057 - 0.0061651*T-0.00000219*T^2 - 0.000000078*T^3; %degs
    planet_coes.w_hat = 48.123691 +0.0291587*T +0.00007051*T^2-0.000000000*T^3; %degs
    planet_coes.L = 304.348665+218.4862002*T+0.00000059*T^2-0.000000002*T^3; %degs
end

%Convert to km:
au = 149597870;
planet_coes.a = planet_coes.a*au;

planet_coes.L = Functions.angle_normalize(planet_coes.L,0); % L
planet_coes.w_hat = Functions.angle_normalize(planet_coes.w_hat,0); % w_hat

planet_coes.h = sqrt(planet_coes.a*uS*(1-planet_coes.ecc^2)); %relative angular momentum
planet_coes.m = Functions.angle_normalize(planet_coes.L - planet_coes.w_hat,0); % mean anomaly in degrees
planet_coes.w = Functions.angle_normalize(planet_coes.w_hat - planet_coes.raan,0); % argument of perihelion
[planet_coes.E,planet_coes.TA] = Functions.ecc_anomaly(planet_coes.m*pi/180,planet_coes.ecc); % eccentric anomaly and true anomaly w/ mean anomaly in rads
planet_coes.E = planet_coes.E*180/pi; % ecc anomaly into degrees
end