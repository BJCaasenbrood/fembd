%% Function Reference
% Click a name to open its reference page (generated from its help text),
% or type |doc fembd.Name| in the Command Window.
%
%% Model
%
% * <matlab:doc('fembd.Body') fembd.Body> - elastic or rigid body: mesh, material, K and M
% * <matlab:doc('fembd.Body.rigidSurface') fembd.Body.rigidSurface> - rigid body from an STL surface (no meshing)
% * <matlab:doc('fembd.Body.pointMass') fembd.Body.pointMass> - rigid body from mass properties only
% * <matlab:doc('fembd.Body/setMass') setMass> - set a rigid body's mass properties (e.g. from CAD)
% * <matlab:doc('fembd.Assembly') fembd.Assembly> - bodies, connections, supports; reduced K, M and u = T*q
% * <matlab:doc('fembd.material') fembd.material> - isotropic linear-elastic material
%
%% Assembly methods
%
% * <matlab:doc('fembd.Assembly/add') add> - add a body
% * <matlab:doc('fembd.Assembly/tie') tie> - bonded tie or elastic interface
% * <matlab:doc('fembd.Assembly/couple') couple> - rigid coupling of nodes to a reference point
% * <matlab:doc('fembd.Assembly/joint') joint> - kinematic joint between reference points
% * <matlab:doc('fembd.Assembly/transmission') transmission> - elastic or rigid ratio between joints or bodies (belts, gears)
% * <matlab:doc('fembd.Assembly/pose') pose> - move a joint to a finite position
% * <matlab:doc('fembd.Assembly/fix') fix> - supports
% * <matlab:doc('fembd.Assembly/dof') dof> - global dof indices
% * <matlab:doc('fembd.Assembly/input') input> - named input (actuator, load)
% * <matlab:doc('fembd.Assembly/output') output> - named output (sensor)
% * <matlab:doc('fembd.Assembly/statespace') statespace> - state-space model with the inputs and outputs
% * <matlab:doc('fembd.Assembly/ioMatrices') ioMatrices> - inputs and outputs in the full space
% * <matlab:doc('fembd.Assembly/qdof') qdof> - positions of dofs in q
% * <matlab:doc('fembd.Assembly/system') system> - reduced K, M and T
% * <matlab:doc('fembd.Assembly/assemble') assemble> - unconstrained K and M
% * <matlab:doc('fembd.Assembly/constraints') constraints> - the map u = T*q
% * <matlab:doc('fembd.Assembly/plot') plot> - plot (deformed) bodies
% * <matlab:doc('fembd.Assembly/animate') animate> - slow-motion mode animation
%
%% Analysis
%
% * <matlab:doc('fembd.modes') fembd.modes> - natural frequencies and mode shapes
% * <matlab:doc('fembd.craigBampton') fembd.craigBampton> - Craig-Bampton and Guyan reduction
%
%% Meshing and I/O
%
% * <matlab:doc('fembd.meshSTL') fembd.meshSTL> - tetrahedral mesh of an STL with Gmsh
% * <matlab:doc('fembd.readMsh') fembd.readMsh> - read an ASCII Gmsh .msh file
% * <matlab:doc('fembd.gmsh') fembd.gmsh> - run Gmsh
% * <matlab:doc('fembd.boxSTL') fembd.boxSTL>, <matlab:doc('fembd.cylinderSTL') fembd.cylinderSTL> - simple test geometry
