% FEMBD  Finite elements and multibody dynamics for STL assemblies
%   Documentation: Help browser, Supplemental Software    Manual: docs/fembd.md
%
% Model
%   Body          - Elastic or rigid body: Gmsh Tet4/Tet10 mesh, material, K and M
%   Assembly      - Bodies, ties, couplings, joints, supports; reduced K, M, u = T*q
%   Model         - Nonlinear multibody model (large rotations): simulate, animate
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
%   steel = fembd.material(210e9, 0.3, 7850);
%   b = fembd.Body('part.stl', steel, Size = 5, Scale = 1e-3);
%   a = fembd.Assembly();  a.add(b);  a.fix(1, @(X) X(:,3) < 1e-9);
%   [K, M, T] = a.system();  [V, f] = fembd.modes(K, M, 6);
