# Tutorials

Runnable MATLAB scripts that explain the toolbox step by step. Each one builds its own geometry,
and most check their results against theory. Run them section by section in the editor
(`Ctrl+Enter`), or render one to HTML with `publish('tut01_firstModel')`. Gmsh must be set up
first (see the main README).

| # | Tutorial | Covers |
|---|---|---|
| 1 | [`tut01_firstModel`](tut01_firstModel.m) | The basic workflow: STL → body → assembly → K, M → modes. Cantilever vs Euler–Bernoulli. |
| 2 | [`tut02_meshing`](tut02_meshing.m) | `meshSTL` options, Tet4 vs Tet10 convergence, curved surfaces, meshes from the Gmsh GUI, `Scale` and `transform`. |
| 3 | [`tut03_materials`](tut03_materials.m) | Isotropic and orthotropic materials, body matrices, mass and inertia from M, consistent vs lumped mass. |
| 4 | [`tut04_ties`](tut04_ties.m) | Non-conforming parts, bonded ties, master and slave, `Tol`, partial bonding with `Select`, elastic interfaces. |
| 5 | [`tut05_loadsStatics`](tut05_loadsStatics.m) | Dof layout, `dof` and `qdof`, supports, point loads, gravity, static solution, reactions. |
| 6 | [`tut06_rigidBodies`](tut06_rigidBodies.m) | Rigid bodies: reference dofs, mass properties, rigid tip mass, supports, `RefPoint`. |
| 7 | [`tut07_joints`](tut07_joints.m) | Couplings on elastic bodies, joint types and axes, joint springs, locking and loading joints, mechanism modes. |
| 8 | [`tut08_dynamics`](tut08_dynamics.m) | Rayleigh damping, frequency response (direct and modal), Newmark time integration, state space. |
| 9 | [`tut09_reduction`](tut09_reduction.m) | Guyan and Craig–Bampton reduction, convergence with the number of modes, recovering the full field, reduced state space. |

The full examples in `../examples` combine these features. The reference manual is
[`../docs/fem.md`](../docs/fem.md).
