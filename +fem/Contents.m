% FEM  Linear-elastic finite elements for STL assemblies
%   Manual: docs/fem.md    Tutorials: tutorials/README.md
%
% Model
%   Body          - Elastic or rigid body: Gmsh Tet4/Tet10 mesh, material, K and M
%   Assembly      - Bodies, ties, couplings, joints, supports; reduced K, M, u = T*q
%   material      - Isotropic linear-elastic material (D, rho)
%
% Analysis
%   modes         - Lowest natural frequencies and mass-normalised modes
%   craigBampton  - Craig-Bampton (or Guyan) reduction
%
% Meshing and I/O
%   meshSTL       - Tetrahedral mesh of a closed STL surface with Gmsh
%   readMsh       - Read tetrahedra from an ASCII Gmsh .msh file (v2.2 / v4.1)
%   gmsh          - Run the Gmsh executable ($GMSH or gmsh on the PATH)
%   boxSTL        - Write an axis-aligned box as STL
%   cylinderSTL   - Write a z-aligned cylinder as STL
%
% Example
%   steel = fem.material(210e9, 0.3, 7850);
%   b = fem.Body('part.stl', steel, Size = 5, Scale = 1e-3);
%   a = fem.Assembly();  a.add(b);  a.fix(1, @(X) X(:,3) < 1e-9);
%   [K, M, T] = a.system();  [V, f] = fem.modes(K, M, 6);
