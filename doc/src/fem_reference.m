%% Function Reference
% Click a name to open its reference page (generated from its help text),
% or type |doc fmbd.Name| in the Command Window.
%
%% Bodies
%
% * <matlab:doc('fmbd.Body') fmbd.Body> - rigid or elastic body: geometry, material, mass properties
% * <matlab:doc('fmbd.Body.elastic') fmbd.Body.elastic> - elastic body meshed with Gmsh
% * <matlab:doc('fmbd.Body.rigid') fmbd.Body.rigid> - rigid body from an STL surface (no meshing)
% * <matlab:doc('fmbd.Body.pointMass') fmbd.Body.pointMass> - rigid body from mass properties only
% * <matlab:doc('fmbd.Body/setMass') setMass> - set a rigid body's mass properties (e.g. from CAD)
% * <matlab:doc('fmbd.Body/massProperties') massProperties> - mass, centre of mass and inertia
% * <matlab:doc('fmbd.material') fmbd.material> - isotropic linear-elastic material
%
%% System: model description
%
% * <matlab:doc('fmbd.System') fmbd.System> - the multibody system
% * <matlab:doc('fmbd.System/add') add> - add a body
% * <matlab:doc('fmbd.System/couple') couple> - rigid coupling of nodes to a reference point
% * <matlab:doc('fmbd.System/tie') tie> - bonded tie or elastic interface between surfaces
% * <matlab:doc('fmbd.System/joint') joint> - kinematic joint, with springs, dampers, Position, Velocity
% * <matlab:doc('fmbd.System/transmission') transmission> - belt, gear or cable between coordinates
% * <matlab:doc('fmbd.System/fix') fix> - supports of reference points, joints and nodes
% * <matlab:doc('fmbd.System/input') input> - named input (actuator, load)
% * <matlab:doc('fmbd.System/output') output> - named output (sensor)
% * <matlab:doc('fmbd.System/update') update> - rebuild after changing a body directly
%
%% System: analysis
%
% * <matlab:doc('fmbd.System/simulate') simulate> - time integration (generalised-alpha)
% * <matlab:doc('fmbd.System/equilibrium') equilibrium> - static equilibrium
% * <matlab:doc('fmbd.System/linearize') linearize> - linear model about a state
% * <matlab:doc('fmbd.System/statespace') statespace> - modal state-space model
% * <matlab:doc('fmbd.System/modes') modes> - natural frequencies and mode shapes
% * <matlab:doc('fmbd.System/energy') energy> - kinetic, gravity, strain and spring energies
% * <matlab:doc('fmbd.System/plot') plot>, <matlab:doc('fmbd.System/animate') animate> - draw a state, play a simulation
% * <matlab:doc('fmbd.System/plotMode') plotMode>, <matlab:doc('fmbd.System/animateMode') animateMode> - draw or animate a mode shape
% * <matlab:doc('fmbd.System/modeShapes') modeShapes> - nodal mode shapes of an elastic body
%
%% Utilities
%
% * <matlab:doc('fmbd.meshSTL') fmbd.meshSTL> - tetrahedral mesh of an STL with Gmsh
% * <matlab:doc('fmbd.readMsh') fmbd.readMsh> - read an ASCII Gmsh .msh file
% * <matlab:doc('fmbd.gmsh') fmbd.gmsh> - run Gmsh
% * <matlab:doc('fmbd.modes') fmbd.modes>, <matlab:doc('fmbd.craigBampton') fmbd.craigBampton> - modes and reduction of K, M
% * <matlab:doc('fmbd.boxSTL') fmbd.boxSTL>, <matlab:doc('fmbd.cylinderSTL') fmbd.cylinderSTL> - simple test geometry
