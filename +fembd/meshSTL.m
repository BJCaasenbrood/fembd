function [X, T] = meshSTL(file, opts)
%MESHSTL Tetrahedral volume mesh of a closed STL surface using Gmsh.
%   [X, T] = fembd.meshSTL(file) returns nodes X (n x 3) and quadratic
%   tetrahedra T (ne x 10, Gmsh ordering).
%
%   Options (name-value):
%     Size          maximum element size in STL units (default: Gmsh default)
%     Order         1 (Tet4) or 2 (Tet10, default)
%     Remesh        true: classify the STL into patches and remesh the
%                   surface (default); false: keep the STL triangles as-is
%     FeatureAngle  angle (deg) used to split the STL into patches (40)
%
%   The Gmsh executable is taken from $GMSH, otherwise 'gmsh' (see fembd.gmsh).

arguments
    file (1,:) char
    opts.Size (1,1) double {mustBePositive} = inf
    opts.Order (1,1) double {mustBeMember(opts.Order, [1 2])} = 2
    opts.Remesh (1,1) logical = true
    opts.FeatureAngle (1,1) double = 40
end

file = strrep(char(java.io.File(file).getCanonicalPath()), '\', '/');
assert(isfile(file), 'fembd:meshSTL', 'File not found: %s', file);

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
    [status, out] = fembd.gmsh(args);
    assert(status == 0 && isfile([base '.msh']), 'fembd:meshSTL', ...
           'Gmsh failed on %s:\n%s', file, out);
    if ~contains(out, 'jac. < 0') || opts.Order == 1, break, end
    if ~straight
        warning('fembd:meshSTL', ['%s: curved Tet10 mesh has inverted elements; ' ...
                'using straight-sided Tet10 instead.'], file);
    end
end

[X, T] = fembd.readMsh([base '.msh']);
end
