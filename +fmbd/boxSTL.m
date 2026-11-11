function boxSTL(file, lo, hi)
%BOXSTL Write an axis-aligned box [lo, hi] (1x3 corners) as a closed STL.
%   Handy for quick test geometry, e.g. fmbd.boxSTL('beam.stl', [0 0 0], [1000 50 50]).

[i, j, k] = ndgrid(0:1);
V = lo + [i(:) j(:) k(:)].*(hi - lo);
F = [1 3 4 2; 5 6 8 7; 1 2 6 5; 3 7 8 4; 1 5 7 3; 2 4 8 6];
stlwrite(triangulation([F(:,[1 2 3]); F(:,[1 3 4])], V), file);
end
