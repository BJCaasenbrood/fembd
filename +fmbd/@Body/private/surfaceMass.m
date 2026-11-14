function [m, c, J, flip] = surfaceMass(X, F, rho)
% Mass, centre of mass and inertia (about c) of a closed triangulated
% surface, summing signed tetrahedra (origin, a, b, c).
o = mean(X, 1);  P = X - o;
a = P(F(:,1),:);
b = P(F(:,2),:);
e = P(F(:,3),:);
s = a + b + e;
d = dot(a, cross(b, e, 2), 2); % 6 x signed volume
V = sum(d)/6;
S1 = sum(d.*s, 1)/24; % int x dV
S2 = (a'*(d.*a) + b'*(d.*b) + e'*(d.*e) + s'*(d.*s))/120; % int x x' dV
flip = V < 0;
if flip, [V, S1, S2] = deal(-V, -S1, -S2); end
m = rho*V;
r = S1/V;
c = o + r;
J = rho*(trace(S2)*eye(3) - S2) - m*(dot(r, r)*eye(3) - r'*r);
end
