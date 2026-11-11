# Tutorials

Runnable MATLAB scripts that explain the toolbox step by step. Each one builds its own geometry,
and most check their results against theory. Run them section by section in the editor
(`Ctrl+Enter`), or render one to HTML with `publish('tut01_firstModel')`. Gmsh must be set up
first (see the main README).

**Modelling and simulation**

| # | Tutorial | Covers |
|---|---|---|
| 1 | [`tut01_firstModel`](tut01_firstModel.m) | The workflow: STL → rigid and elastic bodies → joints, couplings, inputs, outputs → nonlinear simulation with PD control → energy, animation → linear modes and state space. |
| 2 | [`tut02_meshing`](tut02_meshing.m) | `meshSTL` options, Tet4 vs Tet10 convergence, curved surfaces, meshes from the Gmsh GUI, `Scale` and `transform`. |
| 3 | [`tut03_materials`](tut03_materials.m) | Isotropic and orthotropic materials, body matrices, mass and inertia from M, consistent vs lumped mass. |
| 4 | [`tut04_rigidBodies`](tut04_rigidBodies.m) | Rigid bodies: reference dofs, mass properties, rigid tip mass, supports, `RefPoint`, mesh-free rigid bodies. |
| 5 | [`tut05_joints`](tut05_joints.m) | Couplings on elastic bodies, joint types and axes, joint springs, locking and loading joints, mechanism modes, posing, closed loops, belt transmissions. |
| 6 | [`tut06_ties`](tut06_ties.m) | Non-conforming parts, bonded ties, master and slave, `Tol`, partial bonding with `Select`, elastic interfaces. |

**Linear analysis about a pose**

| # | Tutorial | Covers |
|---|---|---|
| 7 | [`tut07_modes`](tut07_modes.m) | `a.system()` → K, M → modes. Cantilever vs Euler–Bernoulli. |
| 8 | [`tut08_loadsStatics`](tut08_loadsStatics.m) | Dof layout, `dof` and `qdof`, supports, point loads, gravity, static solution, reactions. |
| 9 | [`tut09_dynamics`](tut09_dynamics.m) | Rayleigh damping, frequency response (direct and modal), Newmark time integration, state space. |
| 10 | [`tut10_reduction`](tut10_reduction.m) | Guyan and Craig–Bampton reduction, convergence with the number of modes, recovering the full field, reduced state space. |

The full examples in `../examples` combine these features. The reference manual is
[`../docs/fembd.md`](../docs/fembd.md).
