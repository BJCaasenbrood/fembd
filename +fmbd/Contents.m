% FMBD  Flexible multibody dynamics for STL assemblies
%
% Model
%   System        - Multibody system: model, simulation and linear analysis
%   Body          - Rigid or elastic body
%   material      - Isotropic linear-elastic material
%
% Utilities
%   modes         - Natural frequencies and mode shapes of K, M
%   craigBampton  - Craig-Bampton reduction
%   meshSTL       - Mesh an STL surface with Gmsh
%   readMsh       - Read a Gmsh .msh file
%   gmsh          - Run Gmsh
%   boxSTL        - Write a box as STL
%   cylinderSTL   - Write a cylinder as STL
%
% Example
%   s = fmbd.System(Gravity = [0 0 -9.81]);
%   s.add(fmbd.Body.pointMass(1, [0.3 0 0], 1e-3*eye(3), Name = 'link'));
%   s.joint('ground', 'link', 'revolute', Axis = [0 1 0], Name = 'hinge', Position = 0.5);
%   s.output('th', 'hinge');
%   [t, y, x] = s.simulate([0 2]);  s.animate(t, x);
