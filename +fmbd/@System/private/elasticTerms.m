function [Mb, Fb] = elasticTerms(B, R, eta, nu, g)
% Mass matrix and generalised forces (gravity, centrifugal and Coriolis,
% elastic and damping) of an elastic body in its floating frame.
% Velocities nu = [v (world); omega (body); deta].
ne = B.ne;  xi = [1; eta];  S = B.s*xi;  S1 = B.s(:, 2:end);
Quu = reshape(B.Q9*kron(xi, xi), 3, 3);
J = trace(Quu)*eye(3) - Quu;  J = (J + J')/2;
Mwe = reshape(B.CRx*xi, 3, ne);
Mb = [B.m*eye(3), -R*hat(S), R*S1; (-R*hat(S))', J, Mwe; (R*S1)', Mwe', B.Mee];
Mb = (Mb + Mb')/2;
if nargout < 2, return, end
w = nu(4:6);
de = nu(7:end);
gb = R'*g;
Quz = reshape(B.Qaz*kron(de, xi), 3, 3);
W1 = -Quu*hat(w);  W2 = Quz; % u x (w x y) = w (u.y) - y (u.w)
hv = R*(cross(w, cross(w, S)) + 2*cross(w, S1*de));
hw = w*trace(W1) - W1'*w + 2*(w*trace(W2) - W2'*w);
Qju = reshape(B.Qju*xi, 9, ne); % columns: vec of sum_b xi_b Q_{j+1,b}
Qjz = reshape(B.Qjz*de, 9, ne); % columns: vec of sum_k de_k Q_{j+1,k+1}
crz = Qjz([8 3 4], :) - Qjz([6 7 2], :); % cr(Q) = [Q23-Q32; Q31-Q13; Q12-Q21]
he = (kron(w, w)'*Qju)' - (w'*w)*(Qju(1,:) + Qju(5,:) + Qju(9,:))' - 2*(w'*crz)';
Fb = [B.m*g - hv; cross(S, gb) - hw; S1'*gb - he - B.Kc*eta - B.Dc*de];
end
