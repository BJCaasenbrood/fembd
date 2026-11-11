%% fmbd: Flexible Multibody Dynamics
% *fmbd* models assemblies of STL parts as flexible multibody systems.
% Rigid parts need only their surface or mass properties; elastic parts are
% meshed with Gmsh and reduced to a few modes in a moving frame. Parts are
% connected by joints, couplings, ties, springs, dampers and transmissions.
% One nonlinear model gives simulations of large motions (with controllers
% in the loop), equilibria, and linear models about any operating point.
%
%   s = fmbd.System(Gravity = [0 0 -9.81]);
%   s.add(fmbd.Body.rigid('link.stl', 7850, Scale = 1e-3));
%   s.add(fmbd.Body.elastic('tool.stl', fmbd.material(70e9, 0.33, 2700), Scale = 1e-3));
%   s.joint('ground', 'link', 'revolute', Name = 'shoulder', Position = 0.3);
%   s.joint('link', s.couple('tool', @(X) X(:,1) < 0.201), 'fixed');
%   s.input('tau', 'shoulder');  s.output('th', 'shoulder');
%
%   [t, y, x] = s.simulate([0 1], Input = @(t, y) 20*(1 - y.th));
%   s.animate(t, x);
%   x0 = s.equilibrium();  sys = s.statespace(x0, Modes = 20);
%
% Installation and Gmsh setup: see README.md in the toolbox folder. The
% reference manual (model description, analyses, formulation, recipes) is
% <matlab:open(fullfile(fileparts(fileparts(which('fmbd.System'))),'docs','fmbd.md')) docs/fmbd.md>.
%
%% Tutorials
% Runnable scripts in the |tutorials| folder, each checked against theory
% where possible.
%
% # <tut01_firstModel.html Your First Model>: a two-link arm with an
%   elastic forearm: simulation with PD control, then linear analysis
% # <tut02_meshing.html Meshing>: |meshSTL| options, Tet4 vs Tet10, Gmsh
%   GUI meshes, |Scale| and |transform|
% # <tut03_materials.html Materials and Mass>: orthotropic materials, mass
%   properties, Rayleigh damping
% # <tut04_rigidBodies.html Rigid Bodies>: from a mesh, a surface or mass
%   properties; tip masses; |setMass|
% # <tut05_joints.html Couplings, Joints and Transmissions>: joint types,
%   |Position|, large motions, closed loops, belts and gears
% # <tut06_ties.html Ties and Supports>: non-matching meshes, elastic
%   interfaces, node supports
% # <tut07_linear.html Linear Analysis>: equilibrium, modes at several
%   operating points, state space, a controller checked nonlinearly
% # <tut08_elastic.html Elastic Bodies in Motion>: a flexible pendulum,
%   modes, damping, the deformation field
%
%% Examples
%
% * <demoAssembly.html STL Assembly with a Rigid Payload>: elastic parts
%   meshed separately and tied, a rigid plate; modes and a state space
% * <demoSerialArm.html Two-Link Arm with Joints>: elastic links on
%   revolute joints, linearised at several elbow angles
% * <demoParallelFlexure.html Parallel Leaf-Spring Guide>: a flexure stage
%   against beam theory
%
%% Reference
% <fem_reference.html Function Reference>, or |doc fmbd.System|,
% |doc fmbd.Body|.
