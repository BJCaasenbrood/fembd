function [R, t] = jointMotion(J, v)
%JOINTMOTION Finite rigid motion of joint J at coordinates v (its free
%   coordinates: translations in the joint frame, then a rotation vector):
%   x -> R*(x - p) + p + t. Joint coordinates are measured from the
%   configuration in which the joint was created.
s = zeros(6, 1);  s(J.free) = v;
t = J.R*s(1:3);  R = expso3(J.R*s(4:6));
end
