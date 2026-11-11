% FMBD  Flexible multibody dynamics for STL assemblies
%
% Model
%   System        - Bodies, couplings, joints, transmissions, supports, inputs
%                   and outputs; simulate, energy, plot, animate
%   Body          - Rigid or elastic body: Body.rigid, Body.elastic, Body.pointMass
%   material      - Isotropic linear-elastic material (D, rho)
%
% Utilities
%   modes         - Lowest natural frequencies and mass-normalised modes
%   craigBampton  - Craig-Bampton (or Guyan) reduction
%   meshSTL       - Tetrahedral mesh of a closed STL surface with Gmsh
%   readMsh       - Read tetrahedra from an ASCII Gmsh .msh file (v2.2 / v4.1)
%   gmsh          - Run the Gmsh executable ($GMSH or gmsh on the PATH)
%   boxSTL        - Write an axis-aligned box as STL
%   cylinderSTL   - Write a z-aligned cylinder as STL
%
% Example
%   s = fmbd.System(Gravity = [0 0 -9.81]);
%   s.add(fmbd.Body.pointMass(1, [0.3 0 0], 1e-3*eye(3), Name = 'link'));
%   s.joint('ground', 'link', 'revolute', Axis = [0 1 0], Name = 'hinge', Position = 0.5);
%   s.output('th', 'hinge');
%   [t, y, x] = s.simulate([0 2]);  s.animate(t, x);
