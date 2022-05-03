%% ode45 Euler, angular velocity, quaternion propogation function

function [istate] = fiFunc(t,state,I,T)
fi = state(1);
theta = state(2);
psi = state(3);

wx = state(4);
wy = state(5);
wz = state(6);

w = [wx;wy;wz];

Ix = I(1);
Iy = I(2);
Iz = I(3);

Tx = T(1);
Ty = T(2);
Tz = T(3);

epsX = state(7);
epsY = state(8);
epsZ = state(9);
epsilon = [epsX epsY epsZ];

euler = (1/cos(theta)) *     [cos(theta) sin(fi)*sin(theta) cos(fi)*sin(theta);...
                              0          cos(theta)*cos(fi) -sin(fi)*cos(theta);...
                              0          sin(fi)             cos(fi)]                * w;

wxdot = (Tx -(Iz - Iy)*wy*wz)/Ix;
wydot = (Ty - (Ix - Iz)*wx*wz)/Iy;
wzdot = (Tz - (Iy - Ix)*wx*wy)/Iz;

wnew = [wxdot; wydot; wzdot];

epsilonX = [0 -epsZ epsY;...
            epsZ 0 -epsX;...
            -epsY epsX 0];
        
eta = state(10);

epsilonDot = 0.5*(eta*eye(3)+epsilonX)*w;
etaDot = -0.5*epsilon*w;

istate = [euler; wnew; epsilonDot; etaDot];
end