function cylinderSTL(file, c, r, z0, z1, n)
%CYLINDERSTL Write a z-aligned cylinder as a closed STL.
%   fmbd.cylinderSTL(file, c, r, z0, z1) with axis through c = [x y], radius
%   r, from z0 to z1. Optional n: number of facets around (default 64).

if nargin < 6, n = 64; end
th = 2*pi*(0:n-1)'/n;
k = (1:n)';
k1 = mod(k, n) + 1;
ring = c(:)' + r*[cos(th) sin(th)];
V = [ring, z0 + 0*th; ring, z1 + 0*th; c(:)', z0; c(:)', z1];
F = [k k1 n+k1; k n+k1 n+k; 0*k+2*n+1 k1 k; 0*k+2*n+2 n+k n+k1];
stlwrite(triangulation(F, V), file);
end
