# Tutorials

These scripts introduce the toolbox one topic at a time. Each builds its own geometry, and most
compare their results with a hand calculation. Run them section by section in the editor
(Ctrl+Enter) or render them with `publish`. Gmsh has to be set up first (see the main README).

1. [`tut01_firstModel`](tut01_firstModel.m): a two-link arm with an elastic forearm, from STL to
   simulation and linear analysis.
2. [`tut02_meshing`](tut02_meshing.m): meshing options, Tet4 against Tet10, meshes from the Gmsh
   GUI, placing parts.
3. [`tut03_materials`](tut03_materials.m): isotropic and orthotropic materials, body matrices,
   mass properties and damping.
4. [`tut04_rigidBodies`](tut04_rigidBodies.m): rigid bodies from a mesh, a surface or mass
   properties, and CAD mass properties.
5. [`tut05_joints`](tut05_joints.m): joint types and springs, initial positions, a pendulum at
   large amplitude, a four-bar loop, belts and gears.
6. [`tut06_ties`](tut06_ties.m): ties between non-matching meshes, compliant interfaces, node
   supports.
7. [`tut07_linear`](tut07_linear.m): equilibrium, modes over several poses, state-space models and
   a controller check.
8. [`tut08_elastic`](tut08_elastic.m): a flexible pendulum, the number of modes, damping and the
   deformation field.

The scripts in `../examples` combine these features in larger models.
