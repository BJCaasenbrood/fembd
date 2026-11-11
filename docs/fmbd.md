# fmbd — reference manual

`fmbd` models assemblies of STL parts as flexible multibody systems. Rigid parts need only their
surface or mass properties; elastic parts are meshed with Gmsh (Tet4/Tet10) and reduced to a few
modes in a moving frame. Parts are connected by joints, couplings, ties, springs, dampers and
transmissions. One nonlinear model gives:

- **simulation** of large motions with exact joint kinematics, gravity and controllers in the loop;
- **equilibria** under gravity, springs and constant inputs;
- **linear models** about any operating point: modes, state space for control design, frequency
  responses. The stiffness includes gravity, constraint forces and preload, so the models are
  valid at every pose of the mechanism.

New to the toolbox? Start with the tutorials in [`../tutorials`](../tutorials/README.md).

Contents: [1 Workflow](#1-workflow) · [2 Model description](#2-model-description) ·
[3 Analysis](#3-analysis) · [4 Formulation](#4-formulation) · [5 Recipes](#5-recipes) ·
[6 Limitations](#6-limitations)

---

## 1 Workflow

```matlab
alu  = fmbd.material(70e9, 0.33, 2700);                            % E, nu, rho (SI)
link = fmbd.Body.rigid('link.stl', 7850, Scale = 1e-3);            % rigid, STL in mm
tool = fmbd.Body.elastic('tool.stl', alu, Size = 3, Scale = 1e-3, Modes = 10);

s = fmbd.System(Gravity = [0 0 -9.81]);
s.add(link);  s.add(tool);
s.joint('ground', 'link', 'revolute', Point = [0 0 0], Name = 'shoulder', Position = 0.3);
s.joint('link', s.couple('tool', @(X) X(:,1) < 0.201), 'fixed');   % bolt the tool's flange
s.input('tau', 'shoulder');
s.output('th', 'shoulder');  s.output('w', 'shoulder', Type = 'velocity');

% simulation: a PD-controlled move
[t, y, x] = s.simulate([0 1], Input = @(t, y) 80*(1 - y.th) - 8*y.w);
s.animate(t, x, Scale = 100);                    % deformation magnified

% linear analysis about the equilibrium (gravity sag; the free joint stays at its Position)
x0  = s.equilibrium();
[f, V] = s.modes(x0);                            % Hz; s.plotMode(x0, V(:,k))
sys = s.statespace(x0, Modes = 20, Format = "ss");
```

A model goes through three stages:
1. **Bodies** (`fmbd.Body`): geometry plus a material (elastic) or mass properties (rigid).
2. **System** (`fmbd.System`): bodies, connections, supports, inputs and outputs. The bodies are
   defined in the *reference configuration*, in which all joint coordinates are zero.
3. **Analysis**: the system is built on first use (frames, modes, constraints, the initial
   configuration from the joint `Position`s) and rebuilt after any change made through its methods.
   After changing a body directly (`b.Modes`, `b.Damping`, `b.setMass`), call `s.update()`.

---

## 2 Model description

### Bodies
| Constructor | Body |
|---|---|
| `b = fmbd.Body.elastic(file, mat, ...)` | Elastic, meshed with Gmsh |
| `b = fmbd.Body.rigid(file, rho, ...)` | Rigid, from its closed STL surface (no meshing): exact mass, centre of mass and inertia |
| `b = fmbd.Body.pointMass(m, c, J, ...)` | Rigid, from mass properties only (`J` about `c`); no surface, so attach it with a joint |
| `b = fmbd.Body(src, mat, ...)` | General: `src` is an `.stl`, an ASCII `.msh` (v2.2/v4.1) or a struct with `X`, `T`; `Rigid = true` for a meshed rigid body |

| Option | Default | Meaning |
|---|---|---|
| `Name` | file name | Body name |
| `Scale` | 1 | Coordinate factor, e.g. `1e-3` for mm → m |
| `Size`, `Order`, `Remesh`, `FeatureAngle` | Gmsh default, 2, true, 40 | Meshing (`fmbd.meshSTL`) |
| `Modes` | 10 | Fixed-interface vibration modes of an elastic body |
| `Damping` | `[0 0]` | Rayleigh damping `[α β]`: `D = αM + βK`; `β = 2ζ/ω` gives ratio ζ at ω |
| `RefPoint` | centroid / centre of mass | Reference point of a rigid body |

Methods: `b.matrices()` (K, M of the mesh, also `Lumped`, `MassOnly`), `b.massProperties()`,
`b.setMass(m, CenterOfMass=, Inertia=, InertiaPoint=, Axes=, LengthScale=)` (e.g. CAD values),
`b.transform(R, t)`, `b.nodes(sel)`. Materials: `fmbd.material(E, nu, rho)` gives `D` (Voigt
`[xx yy zz yz xz xy]`, engineering shear) and `rho`; any symmetric positive-definite 6×6 `D` works.

### System
| Method | Description |
|---|---|
| `s = fmbd.System(Gravity=)` | New system (`s.Gravity` can be changed later) |
| `s.add(body)` | Add a body; a rigid body is also a reference point with the body's name |
| `name = s.couple(body, sel, Name=, Point=)` | Rigid coupling (RBE2): the selected nodes of an elastic body follow a new reference point. The first coupling of a body carries its frame |
| `n = s.tie(master, slave, Tol=, Stiffness=, Select=)` | Bond two surfaces (non-matching meshes). Two elastic bodies: rigid tie or elastic interface `Stiffness = [kn kt]` (per area); they become one elastic group. Elastic and rigid: the tied nodes become a coupling with a fixed joint |
| `name = s.joint(A, B, type, Point=, Axis=, Frame=, Free=, Stiffness=, Damping=, Position=, Velocity=, Name=)` | Joint between reference points (`'ground'` for A). Types: `fixed revolute prismatic cylindrical spherical planar free`, or a `Free` mask `[ux uy uz rx ry rz]` in the joint frame (local z = `Axis`) |
| `s.transmission(ids, c, Dirs=, Stiffness=, Damping=, Name=)` | Belt/gear/cable: the stretch `δ = Σ cᵢ sᵢ` of joint coordinates or reference-point coordinates (`Dirs` 1:6, absolute). Elastic with `Stiffness`, otherwise rigid (`δ = 0`) |
| `s.fix(ref, dirs)` | Hold a reference point in `dirs` (1:6, world axes) |
| `s.fix(joint, dirs)` | Hold joint coordinates at their `Position` |
| `s.fix(body, sel, dirs)` | Clamp nodes of an elastic body (`dirs` 1:3); its group is grounded, so the supports must hold it |
| `s.input(name, id, sel, dirs)` | Input signal(s): a force/torque on joint coordinates, a force (1:3) or moment (4:6, world) on a reference point, or a total force on a node set (split over the nodes) |
| `s.output(name, id, sel, dirs, Type=)` | Output signal(s): joint coordinates, reference-point coordinates or the average displacement of a node set; `Type` `displacement`, `velocity` or `acceleration` |

**Signals** are one per direction, named `name` or `name_dir`; `s.InputName`, `s.OutputName`.

**Conventions.**
- In the reference configuration (where the bodies are defined), joint coordinates,
  reference-point coordinates and transmission stretches are zero, and joint springs and elastic
  transmissions are relaxed.
- Joint rotations are rotation vectors in the joint frame (for a revolute joint: the angle);
  translations are in the joint frame on A. Reference-point coordinates are absolute: the
  displacement of the point and the rotation vector of its frame.
- The **initial configuration** has the joint coordinates given by `Position`; the constraints
  and the positions are solved exactly, the other joint coordinates stay as close to zero as the
  constraints allow (closed loops adjust), elastic transmissions start unstretched unless all
  their coordinates are set, and elastic bodies start undeformed. `Velocity` sets the initial
  joint rates.

---

## 3 Analysis

| Method | Description |
|---|---|
| `[t, y, x] = s.simulate(tspan, Step=, Input=, Rho=, Tol=, Initial=)` | Time integration. `Input = @(t, y)` returns the input signals from the outputs of the previous step (a sampled controller). `y`: struct of output signals. `x`: states (`r`, `R`, `e`, `nu`, `u` per step). `Initial`: start from a state (e.g. an equilibrium). `Rho` (0.8): numerical damping of high frequencies, 1 = none |
| `[x0, info] = s.equilibrium(Input=, Initial=, Tol=, MaxIter=)` | Static equilibrium under gravity, springs and constant inputs; free mechanisms stay where they are. The inputs are kept in `x0.u`. Newton with line search; stops at a plateau (accepted at round-off level, otherwise a warning: e.g. a free coordinate with a net force). `info`: `converged`, `iterations`, `force`, `constraint` |
| `lin = s.linearize(x0, Input=)` | `M z'' + C z' + K z = B u`, `y = Cz z + Cv z' + Ca z''` in independent coordinates `z` (`dq = T z`) |
| `sys = s.statespace(x0, Modes=, Damping=, StaticCorrection=, Input=, Format=)` | Modal state space `[η; η̇]`: rigid-body/mechanism modes as double integrators, static correction for truncated modes (exact DC gain), modal `Damping` on top of the physical dampers, acceleration feedthrough. `Format = "ss"` (Control System Toolbox) |
| `[f, V] = s.modes(x0, n, Input=)` | Natural frequencies (Hz; unstable directions negative) and mode shapes |
| `E = s.energy(x)` | Kinetic, gravity, strain, spring and total energy per step |
| `s.plot(x, k, Scale=)`, `s.animate(t, x, Rate=, Scale=, Title=)` | Draw a state, or play a simulation (wall-clock; `Rate = 0.1` is 10× slow motion) |
| `s.plotMode(x0, v, Amplitude=)`, `s.animateMode(x0, v, Frequency=, Amplitude=, Rate=, Duration=)` | Draw or animate a mode shape |
| `[Phi, X0] = s.modeShapes(body)` | Nodal mode shapes of an elastic body: its nodes are at `r + R(X0 + Phi η)` |
| `x = s.perturb(x0, dq)` | A state moved along a tangent direction |

Linear models hold about **states at rest**: equilibria, or the initial configuration. An
operating point that is not an equilibrium (e.g. a joint spring preloaded by a `Position`) is
linearised with its instantaneous constraint forces; hold it with inputs (`equilibrium(Input =
u0)`) to model a servo.

---

## 4 Formulation

### 4.1 Coordinates
The system consists of **units**: each rigid body, and each elastic group (elastic bodies tied
together). A unit has a frame with position `r` (world) and rotation `R`, and an elastic unit adds
modal coordinates `η`. Velocities are `ν = [v (world); ω (body frame); η̇]`. A node of an elastic
unit is at

$$x = r + R\,(X_0 + \Phi\,\eta),$$

with `X_0` the node relative to the frame origin. A rigid body's frame is at its centre of mass,
so its mass matrix is constant, `diag(m I, J)`, with gyroscopic forces `ω × Jω`.

### 4.2 Elements and materials
Elastic bodies use isoparametric Tet4 or Tet10 elements (Gmsh node ordering), small-strain
elasticity with Voigt strain `ε = [εxx εyy εzz γyz γxz γxy]` and `σ = D ε`:

$$K_e = \int_{\Omega_e} B^\top D\, B \,d\Omega, \qquad M_e = \int_{\Omega_e} \rho\, N^\top N \,d\Omega ,$$

with collapsed (Duffy) Gauss rules that are exact on straight-sided elements (Tet4: K 1 point, M
27; Tet10: K 27, M 64). Tet4 locks in bending; use Tet10 (the default) unless the mesh is fine.

### 4.3 Elastic units: reduction
For an elastic group, the bodies' `K`, `M` and Rayleigh `D` are combined and the connections are
eliminated:
- **Couplings (RBE2):** coupled nodes follow their coupling, `u = u_c + θ_c × (x − x_c)`.
- **Rigid ties:** each slave surface node within `Tol` of the master surface follows the closest
  master face point, `u_s = Σ N_j(ξ) u_{m,j}` (Tri3/Tri6 shape functions); use the finer mesh as
  slave.
- **Elastic interfaces:** springs on the gap `g = u_s − Σ N_j u_{m,j}`,
  `K += Gᵀ D_s G`, `D_{s,i} = A_i [k_n n nᵀ + k_t (I − n nᵀ)]` with slave nodal areas `A_i`.
- **Frame:** the first coupling is clamped and carries the frame; the other couplings are 6-dof
  interfaces. A group with node supports is grounded instead: the supports are clamped, all
  couplings are interfaces, and the frame is fixed to ground.

The reduced basis is Craig–Bampton,

$$u \approx V \begin{bmatrix} u_b \\ \eta \end{bmatrix},\qquad
  V = \begin{bmatrix} I & 0 \\ \Psi & \Phi \end{bmatrix},\qquad
  \Psi = -K_{ii}^{-1}K_{ib},\qquad (K_{ii} - \omega^2 M_{ii})\Phi = 0 ,$$

with the interface coordinates `u_b`, the static constraint modes `Ψ` and `Modes` fixed-interface
modes `Φ`, plus the residual static shapes `K_{ii}^{-1} f` of node-set inputs (so their static
response is exact). The modal coordinates `η` collect `u_b` and the fixed-interface amplitudes.

### 4.4 Inertia, elastic and gravity forces
The kinetic energy of an elastic unit is `½ ẋᵀ M ẋ` over its nodes. With `U = [X_0, Φ]` and
`ξ = [1; η]`, all terms follow from the invariants `Q_ab = U_aᵀ M_s U_b` (3×3, `M_s` the scalar
node mass matrix) and the first moments `s_a = 1ᵀ M_s U_a`, computed once:
- mass matrix: `M_vv = m I`, `M_vω = −R ŝ(ξ)`, `M_vη = R S₁`, `M_ωω = tr(Q_uu) I − Q_uu`,
  `M_ωη`, `M_ηη = Φᵀ M Φ`, with `Q_uu = Σ ξ_a ξ_b Q_ab`;
- centrifugal and Coriolis forces from `Lᵀ M (ω × (ω × u) + 2 ω × Φ η̇)`;
- elastic and damping forces `−K_c η − D_c η̇` (`K_c = VᵀKV`, `D_c = VᵀDV`);
- gravity `m g`, `s × Rᵀg`, `S₁ᵀ Rᵀ g`.

The cost per step does not depend on the mesh size.

### 4.5 Joints and constraints
Joint ends are **anchors**: a unit frame, or a coupling that moves and rotates with the modes,
`x_c = r + R(x_{c0} + Ψ_t η)`, `R_c = R exp(Ψ_r η)`. A joint has frames `(a, Q_A)` and `(b, Q_B)` on
its two anchors, coincident in the reference configuration. Its coordinates are the translation
`d = Q_Aᵀ R_Aᵀ (p_B − p_A)` and the rotation vector `φ = log(Q_Aᵀ R_Aᵀ R_B Q_B)`. Locked directions
are exact constraints `g(q) = 0`:
- locked translations: components of `d`;
- locked rotations: dot products of frame axes (`e_k^B · e_i^A = 0` for a revolute joint about k),
  valid at any angle.

Supports of reference points are joints to ground; supported joint coordinates, rigid
transmissions (`δ = 0`) and grounded groups add further rows. Redundant rows (e.g. a planar loop
in 3D, or a fixed body that also has a joint to ground) are detected at the start by
rank-revealing QR and dropped. Springs and dampers act on joint coordinates and transmission
stretches, with forces `−Jᵀ(k s + c ṡ)`.

### 4.6 Time integration
The equations `M(q) ν̇ + Gᵀλ = F(q, ν, t)`, `g(q) = 0` are integrated with generalised-α on the
Lie group R³ × SO(3) of each unit (positions add, rotations update as `R exp(Δ)`), with Newton
iterations. Positions satisfy the constraints at every step; velocities are projected onto
`G ν = 0` (mass-orthogonal). `Rho` (spectral radius at infinity) sets the numerical damping of
high frequencies.

### 4.7 Equilibrium and linearisation
`equilibrium` solves `F(q, 0, u₀) = Gᵀλ`, `g(q) = 0` by Newton with minimum-norm steps and a
backtracking line search. It converges when the force residual (relative to the forces) and the
constraint residual (relative to the model size) are below `Tol`; when the residual stops
decreasing, it accepts a round-off-level residual and otherwise stops with a warning.
`linearize` differentiates the same forces: `K = −∂(F − Gᵀλ)/∂q` (central differences on the
tangent space, with the dynamic multipliers λ), `C = −∂F/∂ν` (exact: the forces are at most
quadratic in ν), `B = ∂F/∂u`, and eliminates the constraints with an orthonormal null-space basis
of `G`. `statespace` diagonalises `(K, M)` (symmetric part), keeps `Modes` modes plus static
correction shapes, and adds `2ζω` and the projected `C`.

---

## 5 Recipes

**Operating points.** Linear models at several poses (gain scheduling, robustness checks):
```matlab
for k = 1:numel(th)
    s = buildArm(th(k));                          % e.g. elbow Position = th(k)
    x0 = s.equilibrium(Input = holdingTorques(th(k)));
    S{k} = s.statespace(x0, Modes = 10, Format = "ss");
end
```

**Gravity sag and natural frequencies under load:** `x0 = s.equilibrium(); s.modes(x0)`.

**Frequency response without toolboxes:**
```matlab
sys = s.statespace(x0);
H = @(w) sys.C*((1i*w*eye(size(sys.A)) - sys.A)\sys.B) + sys.D;
```

**Transfer functions and control design** (Control System Toolbox). `zpk`, `tf`, `bode` and
`feedback` need an `ss` object (`Format = "ss"`):
```matlab
sys = s.statespace(x0, Modes = 20, Format = "ss");
p   = tf('s');
G11 = zpk(sys('enc1', 'tau1'));                  % one channel, by signal names
P = sys({'enc1', 'enc2'}, {'tau1', 'tau2'});      % 2x2 plant
C = tf(zeros(2));
for k = 1:2
    LP     = wlp(k)^2/(p^2 + 2*zlp(k)*wlp(k)*p + wlp(k)^2);
    C(k,k) = (KP(k) + KI(k)/p + KD(k)*p)*LP;      % PID with a 2nd-order low-pass
end
CL = feedback(P*C, eye(2));                      % MIMO closed loop
```
- Keep `ss` for the full model: a `tf` of 40+ states is badly conditioned; use `zpk` per channel,
  or reduce first (`balred`).
- Torque to its own joint or motor angle is collocated and robust; feeding back a coordinate
  behind an elastic transmission adds the transmission mode's phase drop.
- Check the design in the nonlinear model: `s.simulate(..., Input = controller, Initial = x0)`.

**Motors with belts.** A rotor (point mass with inertia `Jm`) on a revolute joint, and an elastic
transmission to the driven coordinate; a base-mounted motor driving the second link of an arm
acts on that link's absolute rotation:
```matlab
s.joint('base', 'rotor2', 'revolute', Point = pm, Damping = bm, Name = 'motor2');
s.transmission({'motor2', 'lower'}, [1/N -1], Dirs = {1, 6}, Stiffness = kb, Damping = cb);
```

**Inputs and outputs on node sets:** `s.input('F', 'beam', @(X) X(:,1) > 0.99, 3)` spreads a
force over the nodes (and adds its static shape to the modes); `s.output('a', 'beam', sel, 3,
Type = 'acceleration')` is the average acceleration of the nodes.

---

## 6 Limitations
- The deformation is linear in each unit's frame: no geometric stiffening by centrifugal or
  gravity preload in the parts (it matters for fast-spinning or highly preloaded slender parts),
  no large deformation, no contact.
- An elastic group needs a coupling (its frame) or node supports that hold it.
- Ties between elastic and rigid bodies are rigid (no elastic interface to a rigid body); chained
  ties (a master node that is itself tied) are not supported.
- Linear models are about states at rest.
- Coupled nodes are rigid (RBE2) and stiffen the coupled patch; a distributing coupling (RBE3) is
  not available.
- STL surfaces must be closed for meshing and for `Body.rigid`; only ASCII `.msh` files are read.
- Requires MATLAB R2022a or later and Gmsh 4.

**Verification:** `tests/testFmbd.m` checks the model against analytical solutions (double
pendulum, heavy top, free fall, pendulum and spring equilibria, belt drives, a four-bar's position
analysis), against the linear finite-element solution of the full mesh (gravity sag, modes,
frequency response, ties), and for energy conservation; all simulations converge at second order
in the time step.
