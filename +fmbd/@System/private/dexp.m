function T = dexp(p)
%DEXP Right Jacobian of SO(3): omega = T*dp/dt for R = exp(hat(p)), with
%   dR/dt = R*hat(omega) (omega in the rotated frame).
th = norm(p);  P = hat(p);
if th < 1e-6 % ensure accuracy for small angles
    T = eye(3) - P/2 + (P*P)/6;
else
    T = eye(3) - (1 - cos(th))/th^2*P + (th - sin(th))/th^3*(P*P);
end
end
