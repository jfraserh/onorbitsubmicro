function [istate,raanDot,wpDot] = orbProp(t,state,u,R,J2)
%WORK IN PROGRESS
r0 = state(1:3);
v0 = state(4:6);

coes = Functions.COES(r0,v0,u); % get several orbital elements
vr = coes.vr; %radial velocity
rmag = coes.rmag; %radial magnitude
a = coes.a; %semimajor axis

chi = sqrt(u)*t/a; %universal variable

z = chi^2/a;
if z > 0
    C = (1-cos(sqrt(z)))/z;
    S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z)^3);
elseif z == 0
    C = 1/2;
    S = 1/6;
else
    C = (cosh(sqrt(-z)) - 1)/-z;
    S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(z)^3);
end

% Newton's method 
fratio = (((rmag*vr)/sqrt(u))*chi^2*C+(1-rmag/a)*chi^3*S+rmag*chi-sqrt(u)*t)/(rmag*vr/sqrt(u)*chi*(1-chi^2*S/a)+(1-rmag/a)*(chi^2)*C+rmag);
count = 0;
while abs(fratio) > 1e-8
   chi = chi - fratio; 
   z = chi^2/a;
   if count > 1000
      break; 
   end
   if z > 0
    C = (1-cos(sqrt(z)))/z;
    S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z)^3);
   elseif z == 0
    C = 1/2;
    S = 1/6;
   else
    C = (cosh(sqrt(-z)) - 1)/-z;
    S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(z)^3);
   end
   
   fratio = (((rmag*vr)/sqrt(u))*chi^2*C+(1-rmag/a)*chi^3*S+rmag*chi-sqrt(u)*t)/(rmag*vr/sqrt(u)*chi*(1-chi^2*S/a)+(1-rmag/a)*(chi^2)*C+rmag);
   count = count + 1;
end
f = 1 - chi^2/rmag*C;
g = t - 1/sqrt(u)*chi^3*S;

rnew = f*r0 + g*v0; % new position vector
rmagnew = norm(rnew); % new position

fdot = (sqrt(u)/(rmagnew*rmag))*(chi^3*S/a-chi);
gdot = 1 - chi^2/rmagnew*C;

vnew = fdot*r0 + gdot*v0; % new velocity vector

istate = [rnew vnew];
if nargin > 3
   raanDot = -((3*sqrt(u)*J2*R^2)/(2*(1-coes.eccmag^2)^2*a^(7/2)))*cosd(coes.inc);
   wpDot = raanDot*((5/2)*sind(coes.inc)^2-2)/cosd(coes.inc);
end

end

