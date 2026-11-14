%% Tutorial 2 - Meshing with Gmsh
% Elastic bodies are meshed from STL files with Gmsh. Here we look at the
% meshing options, compare linear (Tet4) and quadratic (Tet10) elements on
% a cantilever, read a mesh made in the Gmsh GUI, and place parts with
% |Scale| and |transform|.

addpath(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(tempdir, 'fmbdTutorials');  if ~isfolder(out), mkdir(out); end

%% 1. Checking Gmsh
% The toolbox calls the Gmsh executable. It uses the GMSH environment
% variable if set, otherwise 'gmsh' on the PATH, e.g.
%   setenv('GMSH', 'C:\gmsh-4.13.1-Windows64\gmsh.exe')
fmbd.gmsh('--version');

%% 2. Meshing an STL directly
% |fmbd.meshSTL| returns nodes X (n x 3) and tetrahedra T (ne x 4 or 10).
% With Remesh = true (default), Gmsh splits the STL into smooth patches,
% using FeatureAngle to detect edges, and remeshes them to the requested
% Size. This helps with CAD STLs, which are full of long thin triangles.
% Remesh = false keeps the STL triangles as the boundary of the mesh.
cyl = fullfile(out, 'cylinder.stl');
fmbd.cylinderSTL(cyl, [0 0], 10, 0, 60); % r = 10, height 60 (mm)
[X1, T1] = fmbd.meshSTL(cyl, Size = 4);
[X2, T2] = fmbd.meshSTL(cyl, Remesh = false);
fprintf('remeshed: %d tets, STL kept as boundary: %d tets\n', size(T1, 1), size(T2, 1));
figure('Color', 'w');
subplot(1, 2, 1);
tetramesh(T1(:,1:4), X1, 'FaceAlpha', 0.1);
axis equal;
title('Remesh = true');
subplot(1, 2, 2);
tetramesh(T2(:,1:4), X2, 'FaceAlpha', 0.1);
axis equal;
title('Remesh = false');

%% 3. Tet4 vs Tet10: a convergence study
% Linear tetrahedra (Tet4) are too stiff in bending ("locking") and need
% very fine meshes. Quadratic Tet10 elements are accurate on coarse meshes
% and are the default. We compute the first cantilever frequency on
% increasingly fine meshes and compare with Euler-Bernoulli. The beam is
% clamped on the nodes of its end face (|fix(body, sel)|).
beam = fullfile(out, 'beam.stl');
fmbd.boxSTL(beam, [0 0 0], [1000 50 50]);
steel = fmbd.material(210e9, 0.3, 7850);
fEB = 1.875104^2/(2*pi)*sqrt(210e9*0.05^4/12/(7850*0.05^2));
sizes = [50 25 12.5];
err = zeros(2, numel(sizes));  ndof = err;
for order = 1:2
    for k = 1:numel(sizes)
        b = fmbd.Body.elastic(beam, steel, Size = sizes(k), Order = order, Scale = 1e-3, Modes = 1);
        s = fmbd.System();
        s.add(b);
        s.fix(b, @(X) X(:,1) < 1e-6);
        err(order, k) = s.modes([], 1)/fEB - 1;  ndof(order, k) = b.nDof;
    end
end
disp(array2table(100*err, 'VariableNames', compose('h=%g', sizes), ...
     'RowNames', {'Tet4 error %', 'Tet10 error %'}));
figure('Color', 'w');
loglog(ndof', abs(err'), 'o-');  grid on
xlabel('mesh dofs');
ylabel('|f_1/f_{EB} - 1|');
legend('Tet4', 'Tet10');
title('Convergence of the first bending frequency');

%% 4. Curved surfaces and Tet10
% Tet10 edge nodes are placed on the curved surface that Gmsh reconstructs
% from the STL. On coarse meshes of small radii this can fold an element;
% meshSTL then remeshes with straight-sided Tet10 and gives a warning. A
% smaller Size keeps the curved geometry.

%% 5. Meshes from the Gmsh GUI
% For full control (local refinement, physical groups, ...), mesh in Gmsh
% yourself and save as ASCII .msh, version 2.2 or 4.1. Pass the .msh file to
% fmbd.Body, or read it with fmbd.readMsh. Here we create such a file
% through fmbd.gmsh, which takes any Gmsh command-line arguments.
geo = fullfile(out, 'block.geo');
writelines(['SetFactory("OpenCASCADE"); Box(1) = {0, 0, 0, 100, 40, 20};' ...
            'Mesh.MeshSizeMax = 8; Mesh.ElementOrder = 2;'], geo);
msh = fullfile(out, 'block.msh');
fmbd.gmsh(sprintf('"%s" -3 -format msh41 -o "%s" -v 0', geo, msh));
block = fmbd.Body.elastic(msh, steel, Scale = 1e-3);
fprintf('%s.msh: %d nodes, %d elements\n', block.Name, block.nNode, size(block.T, 1));

%% 6. Units and placement
% Scale converts coordinates after meshing, so Size stays in file units.
% |transform(R, t)| moves a body rigidly (x -> R*x + t), e.g. to place
% the same part twice. Bodies are defined in the configuration in which
% joint coordinates are zero (see Tutorial 5).
copy = fmbd.Body.elastic(msh, steel, Scale = 1e-3, Name = 'block2');
Rz = [0 -1 0; 1 0 0; 0 0 1]; % 90 degrees about z
copy.transform(Rz, [0.2 0 0]);
s = fmbd.System();
s.add(block);
s.add(copy);
s.fix(block, @(X) X(:,3) < 1e-6);  s.fix(copy, @(X) X(:,3) < 1e-6);
figure('Color', 'w');
s.plot();
title('Two placed copies of one mesh');
