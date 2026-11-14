function R = expso3(p)
%EXPSO3 Rotation matrix exp(hat(p)) of a rotation vector p (Rodrigues).
th = norm(p);  P = hat(p);
if th < 1e-8 % ensure accuracy for small angles
    R = eye(3) + P + P*P/2;
else
    R = eye(3) + sin(th)/th*P + (1 - cos(th))/th^2*(P*P);
end
end
