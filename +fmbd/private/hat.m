function S = hat(v)
%HAT Skew-symmetric matrix: hat(v)*x = cross(v, x).
S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
