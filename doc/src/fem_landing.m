%% fembd: Linear-Elastic Assemblies
% *fembd* builds linear-elastic finite-element models of assemblies of STL
% parts. Each part is meshed with Gmsh (quadratic Tet10 or linear Tet4
% elements) or treated as rigid. Parts are connected by bonded ties,
% elastic interfaces, couplings and joints. The result is the sparse
% stiffness matrix K and mass matrix M, ready for modal, harmonic,
% transient or reduced-order analysis.
%
%   steel = fembd.material(210e9, 0.3, 7850);
%   arm   = fembd.Body('arm.stl', steel, Size = 4, Scale = 1e-3);
%   a = fembd.Assembly();  a.add(arm);
%   a.fix('arm', @(X) X(:,3) < 1e-9);
%   [K, M, T] = a.system();          % u = T*q
%   [V, f] = fembd.modes(K, M, 10);    % natural frequencies (Hz)
%
%% Get started
%
% * <fem_gettingStarted.html Getting Started>: installation, Gmsh setup and a first model
% * <fem_tutorials.html Tutorials>: nine runnable step-by-step tutorials
% * <fem_examples.html Examples>: complete models
%
%% User guide
%
% * <fem_bodies.html Bodies, Meshing and Materials>
% * <fem_connections.html Connections: Ties, Couplings and Joints>
% * <fem_system.html Dofs, Supports, Loads and the System>
% * <fem_analysis.html Modal Analysis, Dynamics and Reduction>
% * <fem_theory.html Formulation>
%
%% Reference
%
% * <fem_reference.html Function Reference>
% * <matlab:doc('fembd.Body') fembd.Body> and <matlab:doc('fembd.Assembly') fembd.Assembly>
%
%% Features
%
% * STL import and Gmsh meshing (Tet4/Tet10), ASCII .msh 2.2/4.1 import
% * isotropic and anisotropic materials; consistent or HRZ-lumped mass
% * elastic and rigid bodies; rigid bodies from STL surfaces or mass
%   properties only
% * bonded ties and elastic interfaces between non-matching meshes
% * rigid couplings and kinematic joints (revolute, prismatic, spherical,
%   ...) with springs, closed loops, and finite posing of joints
% * exact constraint elimination u = T*q; modes of free systems;
%   Craig-Bampton and Guyan reduction
