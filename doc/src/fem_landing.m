%% fembd: Finite Elements and Multibody Dynamics
% *fembd* models assemblies of STL parts as multibody systems with rigid
% and elastic bodies. Elastic parts are meshed with Gmsh; parts are
% connected by joints, couplings, ties, springs, dampers and
% transmissions. The same model gives a nonlinear simulation for large
% motions (with controllers in the loop) and a linear finite-element
% analysis about a pose (modes, frequency responses, state space,
% Craig-Bampton reduction).
%
%   a = fembd.Assembly();
%   a.add(fembd.Body.rigidSurface('link.stl', 7850, Scale = 1e-3));
%   a.add(fembd.Body('tool.stl', fembd.material(70e9, 0.33, 2700), Scale = 1e-3));
%   a.joint('ground', 'link', 'revolute', Name = 'shoulder');
%   a.joint('link', a.couple('tool', @(X) X(:,1) < 0.201), 'fixed');
%   a.input('tau', 'shoulder');  a.output('th', 'shoulder');
%
%   m = a.nonlinear(Gravity = [0 0 -9.81]);      % nonlinear simulation
%   [t, y, x] = m.simulate([0 1], Input = @(t, y) 20*(1 - y.th));
%   m.animate(t, x);
%
%   sys = a.statespace(Modes = 20);              % linear, about the pose
%
% Installation and Gmsh setup: see README.md in the toolbox folder. The
% reference manual (formulation, API, recipes) is
% <matlab:open(fullfile(fileparts(which('fembd.Assembly')),'..','docs','fembd.md')) docs/fembd.md>.
%
%% Tutorials
% Runnable scripts in the |tutorials| folder, each checked against theory
% where possible.
%
% Modelling and simulation:
%
% # <tut01_firstModel.html Your First Model>: a two-link arm with an
%   elastic forearm, simulated with PD control, then analysed linearly
% # <tut02_meshing.html Meshing>: |meshSTL| options, Tet4 vs Tet10,
%   Gmsh GUI meshes, |Scale| and |transform|
% # <tut03_materials.html Materials and Mass>: orthotropic materials, mass
%   properties, consistent vs lumped mass
% # <tut04_rigidBodies.html Rigid Bodies>: reference dofs, tip masses,
%   |RefPoint|, mesh-free rigid bodies
% # <tut05_joints.html Couplings and Joints>: joint types, springs,
%   locking, posing, closed loops, belt transmissions
% # <tut06_ties.html Ties and Interfaces>: non-matching meshes, bonded
%   ties, partial bonding, elastic interfaces
%
% Linear analysis about a pose:
%
% # <tut07_modes.html Modes>: a cantilever against Euler-Bernoulli
% # <tut08_loadsStatics.html Loads and Statics>: dofs, supports, loads,
%   gravity, reactions
% # <tut09_dynamics.html Dynamics>: damping, frequency response, Newmark,
%   state space
% # <tut10_reduction.html Model Reduction>: Guyan and Craig-Bampton
%
%% Examples
%
% * <demoAssembly.html STL Assembly with a Rigid Payload>: elastic parts
%   meshed separately, a rigid plate on an elastic joint
% * <demoSerialArm.html Two-Link Arm with Joints>: elastic links on
%   revolute joints, posed at several elbow angles
% * <demoParallelFlexure.html Parallel Leaf-Spring Guide>: a flexure stage
%   against beam theory
%
%% Reference
% <fem_reference.html Function Reference>, or |doc fembd.Assembly|,
% |doc fembd.Body|, |doc fembd.Model|.
