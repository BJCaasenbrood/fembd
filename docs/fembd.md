# fembd — reference manual

`fembd` builds linear-elastic finite-element models of STL assemblies in MATLAB. Each part is
meshed separately with Gmsh, as Tet4 or Tet10 elements. Parts can be elastic or rigid and are
connected by bonded ties or elastic joints. The output is the sparse stiffness K and mass M,
for modal, harmonic, transient or reduced-order analysis.

New to the toolbox? Start with the tutorials in [`../tutorials`](../tutorials/README.md).

Contents: [1 Workflow](#1-workflow) · [2 Formulation](#2-formulation) ·
[3 API](#3-api) · [4 Recipes](#4-recipes) · [5 Nonlinear model](#5-nonlinear-model-large-rotations) · [6 Limitations](#6-limitations)

---

## 1 Workflow

```matlab
steel = fembd.material(210e9, 0.30, 7850);                     % E, nu, rho (SI)
base  = fembd.Body('base.stl',  steel, Size = 8, Scale = 1e-3); % mesh in mm, model in m
arm   = fembd.Body('arm.stl',   steel, Size = 4, Scale = 1e-3);
tool  = fembd.Body('tool.stl',  struct('rho', 2700), Scale = 1e-3, Rigid = true);

a = fembd.Assembly();
a.add(base);  a.add(arm);  a.add(tool);
a.tie('base', 'arm');                         % bonded interface
a.tie('arm', 'tool', Stiffness = 1e13);       % elastic joint (N/m^3)
a.fix('base', @(X) X(:,3) < 1e-9);            % clamp the bottom face

[K, M, T] = a.system();                       % reduced system, u = T*q
[V, f] = fembd.modes(K, M, 10);                 % f in Hz
a.plot(T*V(:,1), 0.01);                       % deformed shape of mode 1
```

A model goes through three stages:
1. **Bodies:** mesh plus material, each with its own local K and M.
2. **Assembly:** the bodies with their connections and supports. The *full* dof vector `u` stacks every node dof of every body.
3. **System:** constraints are eliminated with `u = T*q`, leaving independent dofs `q` and `K_q = T'*K*T`, `M_q = T'*M*T`.

---

## 2 Formulation

### 2.1 Elements
Bodies are discretised with isoparametric Tet4 or Tet10 elements, using Gmsh node ordering.
Small-strain elasticity, with Voigt strain `ε = [εxx εyy εzz γyz γxz γxy]` (engineering shear)
and `σ = D ε`:

$$K_e = \int_{\Omega_e} B^\top D\, B \,d\Omega, \qquad
  M_e = \int_{\Omega_e} \rho\, N^\top N \,d\Omega .$$

Integrals use a collapsed (Duffy) Gauss rule that is exact for the polynomial degree of the
integrand on straight-sided elements:
- Tet4: K with 1 point, M with 27 points.
- Tet10: K with 27 points, M with 64 points.

Element matrices are computed in vectorised blocks of elements (`pagemtimes`). Elements are
re-oriented so that det J > 0, and `Body.matrices` warns if a curved element has det J ≤ 0 at a
quadrature point.

**Isotropic material:** `D = blkdiag(λ·1·1ᵀ + 2μ I₃, μ I₃)`. Any symmetric positive-definite
6×6 `D` (orthotropic, anisotropic) can be passed in the material struct instead.

**Lumped mass** (`Lumped = true`) uses HRZ lumping: the diagonal of each consistent element mass
is scaled so the element's total mass is preserved. Plain row-sum lumping is not used because it
gives negative corner masses for Tet10.

### 2.2 Bonded ties (multi-point constraints)
`tie(master, slave)` finds the slave surface nodes that lie within `Tol` of the master surface.
Each such node *s* is projected onto its closest master face, at barycentric point ξ, and
constrained to follow the face:

$$u_s = \sum_j N_j(\xi)\, u_{m,j},$$

where `N` are the Tri3 or Tri6 shape functions of the master face. The constraint is exact, with
no penalty parameter, and lets bodies with non-matching meshes be joined. Make the finer mesh the
slave. When a node is a slave in one tie and a master in another, the chain is resolved
automatically.

### 2.3 Elastic interfaces
`tie(..., Stiffness = [kn kt])` uses the same projection to define the gap `g = G u`, i.e.
`u_s − Σ N_j u_m,j`, and adds a distributed spring layer to K:

$$K \mathrel{+}= G^\top D_s G, \qquad
  D_{s,i} = A_i\left[k_n\, n n^\top + k_t (I - n n^\top)\right].$$

- `kn` and `kt` are stiffnesses per unit area (N/m³), normal and tangential to the master surface.
- `A_i` is the slave nodal area, using HRZ weights for Tri6 faces.
- A scalar `Stiffness` sets `kn = kt`.

### 2.4 Rigid bodies
A rigid body has 6 reference dofs, `[ux uy uz rx ry rz]`, at its reference point `x_r`. Every one
of its nodes follows that motion (small rotations):

$$u_j = u_r + \theta \times (x_j - x_r).$$

These node dofs are eliminated through T, just like tie slaves. Its stiffness is zero.
`T'*M*T` of its mesh mass is its exact 6×6 mass and inertia matrix about `x_r`. In a bonded tie,
the rigid body is always made the master.

A rigid body doesn't need a volume mesh:
- `fembd.Body.rigidSurface(file, rho)` reads the STL triangles and uses them for ties and plotting. Mass, centre of mass `c` and inertia `J_c` are computed exactly from the closed surface (divergence theorem).
- `fembd.Body.pointMass(m, c, J)` has only mass properties, e.g. taken from CAD. It has no surface, so connect it with a joint.

For both, the mass is added to the reference dofs as

$$M_6 = \begin{bmatrix} m I & -m[d]_\times \\ m[d]_\times & J_c - m[d]_\times^2 \end{bmatrix},
  \qquad d = c - x_r .$$

A **coupling** (`a.couple`, like Nastran's RBE2) applies the same relation to a selected node set
of an *elastic* body, such as a bore or an end face. This gives the set a 6-dof reference point at
its centroid, or at `Point`. The coupled patch moves rigidly; the rest of the body stays elastic.
Rigid bodies and couplings are both called *reference points*.

### 2.5 Joints
A joint connects reference point B to reference point A, or to ground, at point `p` with frame
`R = [e_x e_y e_z]`. The motion of a reference point at `p` is
`w(p) = P q`, with `P = [I −[p − x_r]×; 0 I]`. The joint requires

$$P_B\, q_B = P_A\, q_A + \begin{bmatrix} R & 0 \\ 0 & R \end{bmatrix} E\, s
  \quad\Longrightarrow\quad q_B = P_B^{-1}\left(P_A q_A + \hat R E s\right),$$

where `E` selects the released directions (the `Free` mask in the joint frame) and `s` are the
joint coordinates.
- The joint is exact: B's 6 dofs become dependent and `s` becomes independent dofs.
- A `Stiffness` on the released directions adds `diag(k)` on `s`: a torsion spring for a revolute joint, an axial spring for a prismatic one, or a 6-dof bushing for `'free'`.
- A `Damping` on the released directions adds `diag(c)` to the damping matrix `D`: viscous bearing or motor friction, or a damped bushing.
- A released joint without stiffness is a mechanism: `fembd.modes` returns it with f = 0.
- **Closed loops** (e.g. four-bar linkages) are allowed in the linear model. Joints are processed in the order they were created. A joint becomes a *cut joint* if its two ends are already connected through earlier joints, if its B side already follows another joint, or if its B side is supported. Its released coordinates follow from the relative motion, `s = (R̂E)ᵀ(P_B q_B − P_A q_A)`, and its constrained directions give rows `C q = 0` (§2.6). Define driving joints first: earlier joints keep their coordinates independent where possible.
- **Posing.** `a.pose(joint, value)` sets the joint to a finite position: an angle in rad or an offset in m, relative to the pose the joint was created in. Everything on its B side moves rigidly: bodies, couplings, downstream joints and tied bodies. The next `a.system()` then linearises about the new pose. Sweeping a pose shows how the dynamics depend on the configuration. If the moving side is attached elsewhere too (a closed loop or a support), `pose` gives an error: posing closed loops would need a loop-closure solve, which is not implemented yet.

| Type | Free in joint frame |
|---|---|
| `fixed` | — |
| `revolute` | `rz` |
| `prismatic` | `uz` |
| `cylindrical` | `uz rz` |
| `spherical` | `rx ry rz` |
| `planar` | `ux uy rz` |
| `free` | all six (use `Stiffness`) |

### 2.5b Transmissions
`a.transmission(joints, c)` relates coordinates of different joints: a belt, gear, cable or rack. Their
"stretch" is `δ = Σ cᵢ sᵢ`.
- **Elastic** (`Stiffness = k`): a spring on `δ`, adding `K += k g gᵀ` with `g` the row of coefficients. `Damping = c` adds a damper on `δ`: `D += c g gᵀ`.
- **Rigid** (no stiffness): the exact constraint `δ = 0`, removed by the coordinate partitioning of §2.6. The first joint listed (the driver) keeps its coordinate.

Each entry is a joint (`Dirs`: its joint coordinate) or a reference point, i.e. a rigid body or
coupling (`Dirs`: 1:6 = `[ux uy uz rx ry rz]`, absolute). For example, a motor (revolute) driving the
rz of a spherical joint with ratio `N = θ_motor/θ_joint`, or a base-mounted motor driving the
*absolute* rotation of a link through a belt:
```matlab
a.transmission({'motor', 'shoulder'}, [1/N -1], Stiffness = k, Dirs = {1, 3});
a.transmission({'motor', 'lower'},    [1/N -1], Stiffness = k, Dirs = {1, 6});
```
The belt torque `k δ` acts on each joint in proportion to its coefficient, and the motor inertia is seen as
`N² J_m` at the joint. Shaft loads from belt tension are not modelled.

### 2.6 Constraint elimination
All dependent dofs `S` (tie slaves, the nodes of rigid bodies and couplings, and the B sides of
joints) and the supported dofs `F` are eliminated. The independent dofs are `I = all \ (S ∪ F)` and

$$u = T q, \qquad T(I,:) = \text{identity}, \quad T(S,:) = R(:,I), \quad u_F = 0 .$$

`R` holds the constraint coefficients after chains are resolved. The reduced matrices are
`T'KT` and `T'MT`. `a.qdof(d)` converts global dofs `d` into positions in `q`.

**Closed loops.** Cut joints give a few extra rows `C q = 0` in the coordinates of the tree. They
are removed by *coordinate partitioning*:
1. Normalise the rows of `C`.
2. Choose dependent columns greedily: non-joint dofs first, then the coordinates of joints from the most recent to the earliest. A column is accepted only if it raises the rank. Redundant rows, such as the out-of-plane rows of a planar loop modelled in 3D, are dropped this way.
3. Solve `q_d = −C_d^{+} C_i q_i` by least squares, which is exact for consistent `C`.

This gives `q = T₂ q'` and `T = T_tree T₂`. The result stays sparse and exact, and K and M stay
symmetric positive (semi)definite, so all downstream tools work unchanged.

### 2.7 Craig–Bampton reduction
Partition the dofs of `K_q, M_q` into boundary dofs `b` and interior dofs `i`:

$$q \approx V \begin{bmatrix} u_b \\ \eta \end{bmatrix},\qquad
  V = \begin{bmatrix} I & 0 \\ \Psi & \Phi \end{bmatrix},\qquad
  \Psi = -K_{ii}^{-1}K_{ib},\qquad (K_{ii} - \omega^2 M_{ii})\Phi = 0 .$$

- `Ψ` are the static constraint modes.
- `Φ` are the first `n` fixed-interface modes, mass-normalised.
- `K_cb = VᵀK_qV` is block-diagonal: the boundary stiffness, then `diag(ω²)`.
- `n = 0` gives Guyan reduction.
- `K_ii` must be non-singular: `b`, together with the supports, must prevent rigid motion.

### 2.8 Eigenvalues
`fembd.modes` uses shift-invert Lanczos (`eigs`) with a small negative shift, scaled to the softest
dofs. This also works for free or partly free systems, where K is singular: rigid-body modes come
out with f ≈ 0. Mode shapes are normalised so that `VᵀMV = I`.

---

## 3 API

### Package functions
| Function | Description |
|---|---|
| `mat = fembd.material(E, nu, rho)` | Isotropic material: struct with `E, nu, rho, D`. Custom: `struct('D', D, 'rho', rho)`; a rigid body needs only `rho`. |
| `[X, T] = fembd.meshSTL(file, ...)` | Gmsh volume mesh of a closed STL. Options: `Size` (max element size, STL units), `Order` (1 or 2, default 2), `Remesh` (default true; false keeps the STL triangles), `FeatureAngle` (40°). |
| `[X, T] = fembd.readMsh(file)` | Read tets from an ASCII Gmsh `.msh` file, version 2.2 or 4.1. |
| `[status, out] = fembd.gmsh(args)` | Run Gmsh with an argument string. Uses `$GMSH`, otherwise `gmsh` on the PATH. On Linux, MATLAB's own libraries are hidden from Gmsh. |
| `fembd.boxSTL(file, lo, hi)`, `fembd.cylinderSTL(file, c, r, z0, z1)` | Write simple test geometry as STL |
| `[V, f] = fembd.modes(K, M, n)` | Lowest `n` modes: shapes `V` (mass-normalised) and frequencies `f` (Hz). |
| `[Kr, Mr, V, f] = fembd.craigBampton(K, M, b, n)` | Craig–Bampton reduction. `b`: boundary dofs (indices into K). `n`: number of fixed-interface modes. `f`: their frequencies. |

### `fembd.Body`
`b = fembd.Body(src, mat, Name=Value)`. `src` is an `.stl` file, an ASCII `.msh` file, or
`struct('X', X, 'T', T)`.

| Option | Default | |
|---|---|---|
| `Name` | file name | Name used to refer to the body in an assembly |
| `Scale` | 1 | Coordinate factor, applied after meshing (e.g. `1e-3`: mm → m) |
| `Size`, `Order`, `Remesh`, `FeatureAngle` | | Passed to `fembd.meshSTL` |
| `Rigid` | false | Make the body rigid (§2.4) |
| `RefPoint` | centroid | Reference point of a rigid body (after `Scale`) |
| `Damping` | `[0 0]` | Rayleigh damping `[α β]` of an elastic body: `D = αM + βK`. `β = 2ζ/ω` gives ratio ζ at ω (rising with frequency); `α` also damps the body's rigid motion |

| Member | |
|---|---|
| `X`, `T`, `F` | Nodes (n×3), tets (ne×4 or ne×10), outward boundary faces (Tri3 or Tri6) |
| `Mat`, `Rigid`, `Ref`, `nNode`, `nDof` | Material, rigid flag, reference point, counts |
| `[K, M] = b.matrices(Lumped=false)` | Body matrices (3·nNode square; K = 0 if rigid) |
| `b = fembd.Body.rigidSurface(file, rho, Name=, Scale=, RefPoint=)` | Rigid body from an STL surface, without meshing |
| `b = fembd.Body.pointMass(m, c, J, Name=, RefPoint=)` | Rigid body from mass properties only (`J` about `c`); connect it with a joint |
| `b.setMass(m, CenterOfMass=, Inertia=, InertiaPoint=, Axes=, LengthScale=)` | Override a rigid body's mass properties (e.g. from CAD). Inertia can be given about another point or in principal axes, and in mm; it is validated (symmetric, positive, triangle inequality). A meshed rigid body then ignores its mesh mass |
| `M6 = b.rigidMass()`, `MassProps` | 6×6 mass about `Ref`, and `struct(m, c, J)`, of a mesh-free rigid body |
| `b.transform(R, t)` | Rigid placement `x → R x + t` (also moves `Ref`) |
| `idx = b.nodes(sel)` | Node indices for a selector |

### `fembd.Assembly`
| Method | |
|---|---|
| `a.add(body)` | Append a body |
| `ns = a.tie(master, slave, Tol=, Stiffness=, Select=)` | Bonded tie, or elastic interface if `Stiffness` is given. `Tol` defaults to 0.1 × the median face edge of the finer surface. `Select`: selector restricting the slave nodes. Returns the number of tied nodes. |
| `r = a.couple(body, sel, Name=, Point=)` | Rigid coupling (RBE2) of selected nodes of an elastic body to a new reference point; returns its name |
| `j = a.joint(A, B, type, Point=, Axis=, Frame=, Free=, Stiffness=, Damping=, Name=)` | Joint between reference points (rigid bodies, couplings, or `'ground'` for A); returns its name (§2.5). `Point` defaults to B's reference point, `Axis` to `[0 0 1]`. |
| `t = a.transmission(ids, c, Stiffness=, Dirs=, Name=)` | Elastic or rigid ratio between joint coordinates or reference-point dofs, `δ = Σ cᵢ sᵢ` (§2.5b) |
| `a.pose(joint, value)` | Move a joint to a finite angle/offset (relative to its creation pose); B side moves rigidly (§2.5) |
| `a.fix(id, sel, dirs)` | Support selected body nodes in `dirs ⊂ 1:3`, or a reference point (`dirs ⊂ 1:6`) or joint coordinates with `a.fix(id, [], dirs)`. |
| `d = a.dof(id, sel, dirs)` | Global dof indices of body nodes, of a reference point `[ux uy uz rx ry rz]`, or of a joint's coordinates. |
| `a.input(name, id, sel, dirs)` | Named input (joint torque/force, force/moment on a reference point, or total force on a node set) |
| `a.output(name, id, sel, dirs, Type=)` | Named output: displacement (default), velocity or acceleration; node sets are averaged |
| `sys = a.statespace(Modes=, Damping=, StaticCorrection=, Lumped=, Format=)` | Modal state-space model with the inputs and outputs (Recipes) |
| `[Bu, Cu, inName, outName, type] = a.ioMatrices()` | Inputs and outputs in the full space |
| `j = a.qdof(d)` | Positions of global dofs `d` in `q`, so `u(d) = q(j)`. `d` must be independent (not supported, not a tie slave). |
| `m = a.nonlinear(Gravity=, Modes=)` | Nonlinear multibody model (`fembd.Model`, §5): `m.simulate(tspan, Step=, Input=, Rates=, Rho=, Tol=)`, `m.energy(x)`, `m.animate(t, x, Rate=, Scale=, Title=)` |
| `[K, M, D] = a.assemble(Lumped=false)` | Full matrices, including elastic interfaces; `D` holds body (Rayleigh), joint and transmission dampers |
| `[T, I] = a.constraints()` | `u = T*q`, and the independent global dofs `I` |
| `[K, M, T, D] = a.system(Lumped=false)` | Reduced matrices, `T`, and the reduced damping matrix `D = T'·D·T` |
| `h = a.plot(u, scale)` | Boundary surfaces, deformed by the full vector `u` and coloured by \|u\| |
| `a.animate(u, Frequency=, Title=, Amplitude=, Rate=, Duration=)` | Slow-motion animation of a mode or displacement `u`, clock-driven (speed independent of rendering); the title shows the slow-motion factor |
| `Bodies`, `Links`, `Refs`, `Joints`, `Fixes`, `nDof`, `X` | Bodies, interfaces, reference points, joints, supports, number of full dofs, all nodes |

**References.** A body is referred to by its index, `Name` or handle. A reference point or joint
is referred to by its name: a rigid body's reference point has the body's name. **Selectors** are
a function `@(X) logical` on node coordinates (n×3), or node indices.

**Dof layout** of `u`:
1. node dofs `[ux1 uy1 uz1 ux2 …]` for each body, in the order added;
2. `[ux uy uz rx ry rz]` for each reference point, in the order created;
3. the released coordinates of each joint, in the order created.

---

## 4 Recipes

**Inputs, outputs and a state-space model** (recommended for control):
```matlab
a.input('tau1', 'motor1');                               % torque on a joint coordinate
a.input('Fz', 'gripper', @(X) X(:,1) > 0.98, 3);         % total force split over selected nodes
a.output('enc1', 'motor1');                              % joint angle
a.output('tip', 'tool', [], 1:3, Type = 'acceleration'); % accelerometer (3 signals)
sys = a.statespace(Modes = 20, Damping = 0.02);          % struct: A, B, C, D, names, f, Shapes
```
- **Signals:** each input or output is one signal per direction, named `name` or `name_dir`.
- **States:** modal, `[η; η̇]`. Rigid-body and mechanism modes are double integrators (undamped unless a joint or transmission damper acts on them).
- **Damping:** the modal ratio `Damping` plus the body, joint and transmission dampers, projected exactly onto the retained modes (a full, non-proportional modal damping matrix).
- **Static correction** (default) adds the quasi-static response to the inputs to the modal basis, so the DC gain is exact with few modes.
- **Acceleration outputs** get the feedthrough `D = CΦΦᵀB`.
- **`Format = "ss"`** returns an `ss` object (Control System Toolbox).
- **`a.ioMatrices()`** gives the full-space `Bu`, `Cu` for your own formulations.

**Transfer functions and control design** (Control System Toolbox). `zpk`, `tf`, `bode` and
`feedback` need an `ss` object: the default struct gives "Invalid syntax". Use `Format = "ss"`, or
convert the struct with `ss(sys.A, sys.B, sys.C, sys.D, 'InputName', sys.InputName, 'OutputName', sys.OutputName)`.
```matlab
sys = a.statespace(Modes = 20, Format = "ss");
s   = tf('s');
G11 = zpk(sys('enc1', 'tau1'));                  % one channel, by signal names
G11r = tf(balred(sys('enc1', 'tau1'), 6));       % low-order tf in s

% decentralised PID with a 2nd-order low-pass per loop
P = sys({'enc1', 'enc2'}, {'tau1', 'tau2'});      % 2x2 plant
C = tf(zeros(2));
for k = 1:2
    LP     = wlp(k)^2/(s^2 + 2*zlp(k)*wlp(k)*s + wlp(k)^2);
    C(k,k) = (KP(k) + KI(k)/s + KD(k)*s)*LP;
end
CL = feedback(P*C, eye(2));                      % MIMO closed loop, coupling included
margin(C(1,1)*P(1,1))                            % loop-at-a-time margins
```
- **Keep `ss` for the full model.** A `tf` of 40+ states has badly conditioned polynomial
  coefficients: use `zpk` per channel, or reduce first (`balred`) and then convert to `tf`.
  `ss*tf` products stay `ss`, so the loop gain keeps the full accuracy.
- **Proper controllers.** A pure derivative `KD*s` is improper. Roll it off with a low-pass
  (above) or a derivative filter, `pid(KP, KI, KD, Tf)`.
- **Collocation.** Torque to its own joint or motor coordinate (`tau1 → enc1`) is collocated:
  poles and zeros alternate and the loop is robust. Feeding back a coordinate behind an elastic
  transmission (`tau1 → th1`, link side of a belt) is not collocated: the belt mode puts a
  phase drop of 180° in the loop and limits the bandwidth.
- **Damping.** With `Damping = 0` only the physical dampers (body, joint, transmission) act, and lightly damped modes give
  sharp peaks in the margins. A small modal ratio (0.01–0.02) is usually more realistic.

**Loads and outputs.** Build them in the full space and project with `T`:
```matlab
d  = a.dof('arm', @(X) X(:,1) > 0.199, 3);         % z dofs of the arm tip
fu = zeros(a.nDof, 1);  fu(d) = 1/numel(d);        % 1 N total, spread over the tip
fq = T'*fu;                                        % load in q
C  = zeros(1, a.nDof);  C(d(1)) = 1;  C = C*T;     % output: u_z at one tip node
```
**Gravity** (works for rigid bodies too):
```matlab
[~, Mu] = a.assemble();
g  = zeros(a.nDof, 1);  g(3:3:3*size(a.X, 1)) = -9.81;
fq = T'*(Mu*g);
```
**Static solution:** `u = T*(K\fq)`.

**Frequency response** with Rayleigh damping `Dq = α K + β M`:
```matlab
H = arrayfun(@(w) C*((K - w^2*M + 1i*w*Dq)\fq), 2*pi*fHz);
```
**Modal truncation:** `[V, f] = fembd.modes(K, M, n);` then `Kn = diag((2*pi*f).^2)`,
`Mn = eye(n)`, `fn = V'*fq`, `Cn = C*V`.

**Joints between parts:**
```matlab
r1 = a.couple('arm', @(X) vecnorm(X(:,[1 3]) - c, 2, 2) < 6e-3);   % elastic bore
a.joint('base', r1, 'revolute', Point = c3, Axis = [0 1 0], Stiffness = 200);
a.joint('arm', 'gripper', 'prismatic', Axis = [1 0 0]);           % between rigid bodies
s = a.dof('joint1');                                              % joint coordinate(s)
```
Joint coordinates are ordinary dofs: support them (`a.fix('joint1', [])`) to lock the joint, load
them (a torque on a revolute joint), or use them as Craig–Bampton boundary dofs, which gives an
actuator-level model (see `examples/demoSerialArm.m`).

**Craig–Bampton model** with a rigid flange as interface:
```matlab
ib = a.qdof(a.dof('tool'));                        % 6 reference dofs
[Kc, Mc, Vc] = fembd.craigBampton(K, M, ib, 20);     % u = T*Vc*[u_b; eta]
```
**Stiffness matrix at a rigid body** (e.g. a flexure stage's shuttle), by Guyan condensation onto
its 6 dofs:
```matlab
Ks = fembd.craigBampton(K, M, a.qdof(a.dof('shuttle')), 0);   % 6x6, order [ux uy uz rx ry rz]
Cs = inv(Ks);                                               % compliance at the reference point
```
See `examples/demoParallelFlexure.m`.
**State space** (sparse descriptor form, for large models):
```matlab
n = size(K, 1);  Z = sparse(n, n);
E = [speye(n) Z; Z M];  A = [Z speye(n); -K -Dq];  B = [zeros(n, 1); fq];  Cx = [C zeros(1, n)];
```
For a small reduced model (e.g. CB), use `A = [0 I; -M\K -M\D]` with `ss(A, B, C, 0)`.

**Explicit time integration:** use `a.system(Lumped = true)`. M is diagonal unless bonded ties
spread slave mass onto master nodes; elastic joints keep it diagonal.

**Meshing in Gmsh yourself:** save the mesh as ASCII `.msh` and pass it as `src`, or use
`fembd.readMsh` and pass `struct('X', X, 'T', T)`.

---

## 5 Nonlinear model (large rotations)

`m = a.nonlinear(Gravity = [0 0 -9.81])` turns the same assembly into a nonlinear multibody model
(`fembd.Model`) for large rotations, e.g. a robot arm moving through its workspace. The linear
model (`a.statespace`) stays the tool for frequency-domain analysis about one pose.
```matlab
m = a.nonlinear(Gravity = [0 0 -9.81]);
[t, y, x] = m.simulate([0 2], Step = 1e-3, ...
    Input = @(t, y) [KP1*(r1(t) - y.enc1); KP2*(r2(t) - y.enc2)]);   % controller in the loop
plot(t, y.th1);  m.animate(t, x, Rate = 0.5);                         % half speed
```
- **Coordinates:** each rigid body has a frame at its centre of mass: position `r` (world) and
  rotation `R`, with velocities `[v (world); ω (body)]`. The mass matrix is constant
  (`diag(m I, J)`), with gyroscopic forces `ω × Jω`.
- **Elastic bodies** use a floating frame of reference. The frame follows the body's first
  coupling; the deformation in that frame is described by Craig–Bampton modes (first coupling
  clamped, the other couplings as interfaces, plus `Modes` fixed-interface modes, default 10):
  `x = r + R(X₀ + Φη)`. The mass matrix `M(η)`, the centrifugal and Coriolis forces and gravity
  follow from precomputed inertia invariants `Q_ab = U_aᵀ M U_b`, `U = [X₀, Φ]`, so the cost per
  step does not depend on the mesh size. Stiffness and Rayleigh damping are the reduced
  `Vᵀ K V` and `Vᵀ D V`. Couplings move and rotate with the modes, so joints on them deform with
  the body. Large rotations are exact; the deformation is linear (no stiffening by centrifugal
  or gravity preload, which matters only for fast-spinning slender parts).
- **Joints and supports** are exact nonlinear constraints `g(q) = 0` with Lagrange multipliers.
  Locked translations: the components of the joint offset in the frame on A. Locked rotations:
  dot products of frame axes (a revolute joint keeps its axes aligned at any angle). Redundant
  constraints (e.g. a fixed body that also has a joint to ground) are detected at the start
  (rank-revealing QR) and dropped.
- **Forces:** gravity, joint springs and dampers, elastic transmissions, and the inputs (joint
  forces/torques; forces and moments on rigid bodies, world axes). Rigid transmissions and
  supported joints are constraints.
- **Coordinates of joints** are measured from the configuration in which the joint was created
  (as in `pose`): translations in the joint frame, rotations as the rotation vector of B relative
  to A (for a revolute joint: the angle). Joint springs are relaxed there. Transmissions are
  relaxed (elastic) or held (rigid) at the initial configuration. Reference-point outputs are the
  displacement of the point and the rotation vector of the body, from the initial configuration.
- **Time integration:** generalised-α on the Lie group R³ × SO(3) of each body, with Newton
  iterations. Positions satisfy the constraints at every step; velocities are projected onto
  them. `Rho` (default 0.8) sets the numerical damping of high frequencies (1: none).
- **Inputs** come from `Input = @(t, y)`, evaluated with the outputs of the previous step: a
  controller sampled at the time step. Initial joint velocities: `Rates = struct('joint', v)`.
- **Signals:** joint coordinates, reference points (rigid bodies, couplings) and node sets
  (average displacement, world axes; forces split over the nodes). `m.energy(x)` gives the
  kinetic, gravity, strain and spring energies; `m.animate(t, x, Scale = 100)` magnifies the
  deformation.
- **Limits:** an elastic body needs at least one coupling (its frame); ties and supports on
  body nodes are not supported (support a coupling instead); displacement and velocity outputs.

Verified against the Lagrange equations of a double pendulum, the exact (linear) solution of a
damped belt drive at large angles, a rigid transmission under a constant torque, a heavy top
(energy, Euler equations) and free fall; elastic bodies against the linear model (cantilever
with a tip mass on its second coupling, node-set output), a stiff flexible pendulum against the
rigid one, energy conservation of a soft flexible pendulum, and the kinetic energy from the
invariants against the nodal `½ ẋᵀ M ẋ`. All converge at second order in the time step.

---

## 6 Limitations
- Linear elasticity only: small strains, no contact, no geometric stiffness. Large rotations of rigid and elastic bodies: see the nonlinear model (§5).
- Ties project onto the flat corner triangle of each master face. On curved Tet10 interfaces the error is about the chord height.
- If curved Tet10 edge nodes invert an element (coarse mesh, small radius), `meshSTL` re-meshes with straight-sided Tet10 and warns. A smaller `Size` keeps the curved geometry.
- STL surfaces must be closed and watertight for `Remesh = true`. Only ASCII `.msh` files are read.
- Closed kinematic loops are supported in the linear model, but `pose` cannot move them (no loop-closure solver yet). Near singular poses (e.g. a four-bar at dead centre) the partitioning becomes ill-conditioned, as the mechanism itself does.
- Couplings are rigid (RBE2) and stiffen the coupled patch. A distributing coupling (RBE3) is not available.
- `T` is formed explicitly: memory grows with the number of tie and rigid-body nodes. `craigBampton` forms the dense `V`.
- Requires MATLAB R2022a or later and Gmsh 4.

Verification: `tests/runTests.m` covers 81 checks, including Euler–Bernoulli, rigid-body
invariance, exact rigid inertia, tie and spring consistency, joints against 1-dof oscillators, a four-bar against velocity analysis, and
hinged-beam solutions, couplings, posed joints against a double pendulum, mesh-free rigid bodies, Craig–Bampton accuracy, damping, state space, the MSH readers and the nonlinear model (§5).
