function [coes] = COES(rvec,vvec,u)

coes.rmag = norm(rvec);
coes.vmag = norm(vvec);
coes.vr = (dot(rvec,vvec))/coes.rmag;

% Find angular momentum
coes.hvec = cross(rvec,vvec);
coes.hmag = norm(coes.hvec);

% Find eccentricity
coes.eccvec = 1/u*(cross(vvec,coes.hvec) - u*(rvec/coes.rmag));
coes.eccmag = norm(coes.eccvec);

% Find inclination
coes.inc = acosd(coes.hvec(3)/coes.hmag);

% Find node line
coes.node = cross([0 0 1],coes.hvec);
coes.nodemag = norm(coes.node);

% Find RAAN
if coes.node(2) < 0
    coes.raan = 360 - acosd(coes.node(1)/coes.nodemag);
else
    coes.raan = acosd(coes.node(1)/coes.nodemag);
end

% Find argument of periapse
if coes.eccvec(3) < 0
   coes.w = 360 - acosd(dot(coes.node,coes.eccvec)/(coes.nodemag*coes.eccmag));
else
    coes.w = acosd(dot(coes.node,coes.eccvec)/(coes.nodemag*coes.eccmag));
end

% Find true anomaly
if coes.vr > 0
   coes.TA = acosd(dot(coes.eccvec,rvec)/(coes.eccmag*coes.rmag)); 
else
    coes.TA = 360 - acosd(dot(coes.eccvec,rvec)/(coes.eccmag*coes.rmag)); 
end

% Energy
coes.E = coes.vmag^2/2 - u/coes.rmag;

% Find semimajor axis
if coes.eccmag < 1 && coes.eccmag > 0
    coes.a = (coes.hmag^2/u)*(1/(1-coes.eccmag^2));
elseif coes.eccmag == 0
    coes.a = coes.rmag*2;
elseif coes.eccmag == 1
    
else
    coes.a = -(coes.hmag^2/u)*(1/(coes.eccmag^2-1));
end

% Find period if an ellipse
if coes.eccmag < 1 && coes.eccmag > 0
   coes.T = 2*pi/sqrt(u)*(coes.a)^(3/2); 
   coes.E0 = 2*atan(sqrt((1-coes.eccmag)/(1+coes.eccmag))*(tand(coes.TA/2)));
   coes.M = coes.E0 - coes.eccmag*sin(coes.E0);
   coes.tSinceP = coes.M*coes.hmag^3/u^2/(1-coes.eccmag^2)^(3/2);
   coes.n = coes.M/coes.tSinceP;
end

end