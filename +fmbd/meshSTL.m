function [X, T] = meshSTL(file, opts)
%MESHSTL Tetrahedral mesh of a closed STL surface with Gmsh.
%   [X, T] = fmbd.meshSTL(file) returns nodes X and tetrahedra T (Gmsh order).
%   Options: Size (STL units), Order (1 = Tet4, default; 2 = Tet10),
%   Remesh (remesh the surface, default true), FeatureAngle (default 40).

arguments
    file (1,:) char
    opts.Size (1,1) double {mustBePositive} = inf
    opts.Order (1,1) double {mustBeMember(opts.Order, [1 2])} = 1
    opts.Remesh (1,1) logical = true
    opts.FeatureAngle (1,1) double = 40
end

file = strrep(char(java.io.File(file).getCanonicalPath()), '\', '/');
assert(isfile(file), 'fmbd:meshSTL', 'File not found: %s', file);

if opts.Remesh
    surf = sprintf(['ClassifySurfaces{%g*Pi/180, 1, 1, Pi};\n' ...
                    'CreateGeometry;\n'], opts.FeatureAngle);
else
    surf = '';
end
geo = @(straight) sprintf(['Merge "%s";\n%s' ...
                          'Surface Loop(1) = Surface{:};\n' ...
                          'Volume(1) = {1};\n' ...
                          'Physical Volume(1) = {1};\n' ...
                          'Mesh.SecondOrderLinear = %d;\n'], file, surf, straight);
base = tempname;
cleanup = onCleanup(@() delete([base '.*']));
args = sprintf('"%s.geo" -3 -order %d -format msh22 -o "%s.msh" -v 2', ...
               base, opts.Order, base);
if isfinite(opts.Size), args = sprintf('%s -clmax %g', args, opts.Size); end

% Tet10 edge nodes are first placed on the curved surface; if that inverts
% elements, remesh with straight-sided Tet10
for straight = [0 1]
    writelines(geo(straight), [base '.geo']);
    [status, out] = fmbd.gmsh(args);
    assert(status == 0 && isfile([base '.msh']), 'fmbd:meshSTL', ...
           'Gmsh failed on %s:\n%s', file, out);
    if ~contains(out, 'jac. < 0') || opts.Order == 1, break, end
    if ~straight
        warning('fmbd:meshSTL', ['%s: curved Tet10 mesh has inverted elements; ' ...
                'using straight-sided Tet10 instead.'], file);
    end
end

[X, T] = fmbd.readMsh([base '.msh']);
end
