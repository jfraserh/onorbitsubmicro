function [s] = quatmult(p,q)
%QUATMULT Carries out quaternion multiplication of p and q

[~,col] = size(p);
if col ~= 1
    p = p';
end

[~,col] = size(q);
if col ~= 1
    q = q';
end

np = p(4);
ep = p(1:3);

nq = q(4);
eq = q(1:3);

% Instantiate s
s = zeros([4 1]);

s(1:3,1) = np.*eq+nq.*ep+Functions.vec2cross(ep)*eq;
s(4) = np*nq-ep'*eq;
end

