# Tutorials

Runnable MATLAB scripts that explain the toolbox step by step. Each one builds its own geometry,
and most check their results against theory. Run them section by section in the editor
(`Ctrl+Enter`), or render one to HTML with `publish('tut01_firstModel')`. Gmsh must be set up
first (see the main README).

| # | Tutorial | Covers |
|---|---|---|
| 1 | [`tut01_firstModel`](tut01_firstModel.m) | The workflow: STL → rigid and elastic bodies → couplings, joints, inputs, outputs → simulation with PD control → energy, animation → equilibrium, modes, state space. |
| 2 | [`tut02_meshing`](tut02_meshing.m) | `meshSTL` options, Tet4 vs Tet10 convergence, curved surfaces, meshes from the Gmsh GUI, `Scale` and `transform`. |
| 3 | [`tut03_materials`](tut03_materials.m) | Isotropic and orthotropic materials, body matrices, mass properties, Rayleigh damping. |
| 4 | [`tut04_rigidBodies`](tut04_rigidBodies.m) | Rigid bodies from a mesh, an STL surface or mass properties; a rigid tip mass; supports; `RefPoint`; `setMass` from CAD. |
| 5 | [`tut05_joints`](tut05_joints.m) | Couplings, joint types and springs, loading and locking joints, mechanism modes, `Position` and large motions (pendulum period), a closed four-bar loop, belts and gears. |
| 6 | [`tut06_ties`](tut06_ties.m) | Node supports, bonded ties between non-matching meshes, `Tol`, partial bonding with `Select`, elastic interfaces, ties to rigid bodies. |
| 7 | [`tut07_linear`](tut07_linear.m) | Equilibrium under gravity, modes at several operating points, state space (truncation, static correction, acceleration outputs), a controller designed linearly and checked in the nonlinear model. |
| 8 | [`tut08_elastic`](tut08_elastic.m) | A flexible pendulum: energy exchange, convergence with `Modes`, physical vs numerical damping, the deformation field, limitations. |

The examples in `../examples` combine these features. The reference manual is
[`../docs/fmbd.md`](../docs/fmbd.md).
