%% Function Reference
% Click a name to open its reference page (generated from its help text),
% or type |doc fem.Name| in the Command Window.
%
%% Model
%
% * <matlab:doc('fem.Body') fem.Body> - elastic or rigid body: mesh, material, K and M
% * <matlab:doc('fem.Body.rigidSurface') fem.Body.rigidSurface> - rigid body from an STL surface (no meshing)
% * <matlab:doc('fem.Body.pointMass') fem.Body.pointMass> - rigid body from mass properties only
% * <matlab:doc('fem.Body/setMass') setMass> - set a rigid body's mass properties (e.g. from CAD)
% * <matlab:doc('fem.Assembly') fem.Assembly> - bodies, connections, supports; reduced K, M and u = T*q
% * <matlab:doc('fem.material') fem.material> - isotropic linear-elastic material
%
%% Assembly methods
%
% * <matlab:doc('fem.Assembly/add') add> - add a body
% * <matlab:doc('fem.Assembly/tie') tie> - bonded tie or elastic interface
% * <matlab:doc('fem.Assembly/couple') couple> - rigid coupling of nodes to a reference point
% * <matlab:doc('fem.Assembly/joint') joint> - kinematic joint between reference points
% * <matlab:doc('fem.Assembly/transmission') transmission> - elastic or rigid ratio between joints or bodies (belts, gears)
% * <matlab:doc('fem.Assembly/pose') pose> - move a joint to a finite position
% * <matlab:doc('fem.Assembly/fix') fix> - supports
% * <matlab:doc('fem.Assembly/dof') dof> - global dof indices
% * <matlab:doc('fem.Assembly/input') input> - named input (actuator, load)
% * <matlab:doc('fem.Assembly/output') output> - named output (sensor)
% * <matlab:doc('fem.Assembly/statespace') statespace> - state-space model with the inputs and outputs
% * <matlab:doc('fem.Assembly/ioMatrices') ioMatrices> - inputs and outputs in the full space
% * <matlab:doc('fem.Assembly/qdof') qdof> - positions of dofs in q
% * <matlab:doc('fem.Assembly/system') system> - reduced K, M and T
% * <matlab:doc('fem.Assembly/assemble') assemble> - unconstrained K and M
% * <matlab:doc('fem.Assembly/constraints') constraints> - the map u = T*q
% * <matlab:doc('fem.Assembly/plot') plot> - plot (deformed) bodies
%
%% Analysis
%
% * <matlab:doc('fem.modes') fem.modes> - natural frequencies and mode shapes
% * <matlab:doc('fem.craigBampton') fem.craigBampton> - Craig-Bampton and Guyan reduction
%
%% Meshing and I/O
%
% * <matlab:doc('fem.meshSTL') fem.meshSTL> - tetrahedral mesh of an STL with Gmsh
% * <matlab:doc('fem.readMsh') fem.readMsh> - read an ASCII Gmsh .msh file
% * <matlab:doc('fem.gmsh') fem.gmsh> - run Gmsh
% * <matlab:doc('fem.boxSTL') fem.boxSTL>, <matlab:doc('fem.cylinderSTL') fem.cylinderSTL> - simple test geometry
