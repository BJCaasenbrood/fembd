%% Nonlinear Model (Large Rotations)
% The linear model (|a.system|, |a.statespace|) holds for small motions
% about one pose. For large rotations, e.g. a robot arm moving through its
% workspace, the same assembly gives a nonlinear multibody model:
%
%   m = a.nonlinear(Gravity = [0 0 -9.81]);          % fembd.Model
%   [t, y, x] = m.simulate([0 2], Step = 1e-3, ...
%       Input = @(t, y) [KP1*(r1(t) - y.enc1); KP2*(r2(t) - y.enc2)]);
%   plot(t, y.th1);  m.animate(t, x, Rate = 0.5);    % half speed
%
%% Formulation
% Each rigid body has a frame at its centre of mass, with position r
% (world) and rotation R; its velocities are v (world) and omega (body
% frame). Elastic bodies add modal coordinates eta (see below). The
% equations of motion are
%
% $$M \dot\nu + G^\top \lambda = F(q, \nu, t), \qquad g(q) = 0,$$
%
% with the mass matrix M (constant diag(m I, J) for rigid bodies), the
% forces F (gravity, gyroscopic, elastic, joint springs and dampers,
% elastic transmissions, inputs) and the joint constraints g with
% Jacobian G.
%
% * *Joints and supports* are exact nonlinear constraints: locked
%   translations are components of the joint offset in the frame on A;
%   locked rotations are dot products of frame axes, so a revolute joint
%   keeps its axes aligned at any angle. Redundant constraints are dropped
%   at the start (rank-revealing QR).
% * *Rigid transmissions* and supported joints are constraints; elastic
%   transmissions are springs and dampers.
% * *Time integration*: generalised-alpha on the Lie group R3 x SO(3) of
%   each body, with Newton iterations. Positions satisfy the constraints at
%   every step and velocities are projected onto them. |Rho| (default 0.8)
%   damps high frequencies numerically; 1 gives no numerical damping.
%
%% Elastic bodies
% Elastic bodies use a floating frame of reference. The frame follows the
% body's first coupling; the deformation in that frame is described by
% Craig-Bampton modes (first coupling clamped, the other couplings as
% interfaces, plus |Modes| fixed-interface modes, default 10):
%
% $$x = r + R\,(X_0 + \Phi\,\eta).$$
%
% * The mass matrix M(eta), the centrifugal and Coriolis forces and gravity
%   follow from inertia invariants computed once from the FE mass matrix,
%   so the cost per step does not depend on the mesh size.
% * Stiffness and Rayleigh damping are the reduced V'KV and V'DV.
% * Couplings move and rotate with the modes, so joints on them deform
%   with the body.
% * Large rotations are exact; the deformation is linear (no stiffening by
%   centrifugal or gravity preload, which matters only for fast-spinning
%   slender parts).
%
%   m = a.nonlinear(Gravity = [0 0 -9.81], Modes = 10);
%   [t, y, x] = m.simulate([0 1], Step = 5e-4);
%   E = m.energy(x);                           % kinetic, gravity, strain, springs
%   m.animate(t, x, Scale = 100);              % deformation magnified 100x
%
%% Coordinates and signals
%
% * Joint coordinates are measured from the configuration in which the
%   joint was created, as in |a.pose|: translations in the joint frame,
%   rotations as the rotation vector of B relative to A (for a revolute
%   joint the angle). Joint springs are relaxed there.
% * Transmissions are relaxed (elastic) or held (rigid) at the initial
%   configuration.
% * Reference-point outputs are the displacement of the point and the
%   rotation vector of its frame, from the initial configuration.
%   Node-set signals are the average displacement of the nodes (world
%   axes); node-set forces are split over the nodes.
% * Inputs come from |Input = @(t, y)|, evaluated with the outputs of the
%   previous step (a controller sampled at the time step). Initial joint
%   velocities: |Rates = struct('joint', v)|.
%
%% Current limits
% An elastic body needs at least one coupling (its frame). Ties and
% supports on body nodes are not supported (support a coupling instead).
% Outputs are displacements and velocities.
%
%% Verification
% |tests/runTests.m| checks the model against the Lagrange equations of a
% double pendulum, the exact solution of a damped belt drive at large
% angles (linear for parallel axes), a rigid transmission under constant
% torque, a heavy top (energy and Euler equations) and free fall; elastic
% bodies against the linear model (cantilever with a tip mass, node-set
% output), a stiff flexible pendulum against the rigid one, and energy
% conservation of a soft flexible pendulum. All converge at second order
% in the time step.
