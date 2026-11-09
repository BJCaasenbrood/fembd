%% Dofs, Supports, Loads and the System
%
%% Dof layout
% The full displacement vector u contains, in order:
%
% # the node dofs [ux1 uy1 uz1 ux2 ...] of every body, in the order the
%   bodies were added;
% # 6 dofs [ux uy uz rx ry rz] per reference point (rigid bodies and
%   couplings), in creation order;
% # the released coordinates of every joint, in creation order.
%
% |a.dof| returns global dof indices:
%
%   d = a.dof(body, sel, dirs)    % nodes selected by sel, dirs in 1:3
%   d = a.dof(ref, [], dirs)      % reference point, dirs in 1:6
%   d = a.dof(joint)              % joint coordinates
%
% A selector |sel| is a function of the node coordinates, e.g.
% |@(X) X(:,1) < 1e-9|, or a list of node indices.
%
%% Supports
%
%   a.fix(body, sel, dirs)        % nodes (dirs default 1:3)
%   a.fix(ref, [], dirs)          % rigid body or coupling, dirs in 1:6
%   a.fix(joint, [])              % lock a joint
%
%% The system
%
%   [K, M, T] = a.system(Lumped = false)
%
% All constraints are eliminated with u = T*q: tie slaves, nodes of rigid
% bodies and couplings, B sides of joints, cut joints and supports. The
% returned K = T'*Kfull*T and M = T'*Mfull*T are symmetric and sparse. A
% fourth output D = T'*Dfull*T holds body (Rayleigh), joint and transmission dampers.
% |a.assemble()| gives the unconstrained Kfull and Mfull, and
% |a.constraints()| gives T and the independent dofs.
%
% |a.qdof(d)| converts global dofs into positions in q, so
% u(d) = q(a.qdof(d)).
%
%% Inputs, outputs and state space
% Named inputs (actuators, loads) and outputs (sensors), one signal per
% direction, selected as in |a.dof|:
%
%   a.input('tau1', 'motor1');                          % joint torque
%   a.input('Fz', 'gripper', @(X) X(:,1) > 0.98, 3);    % total force on a node set
%   a.output('enc1', 'motor1');                         % joint angle
%   a.output('tip', 'tool', [], 1:3, Type = 'acceleration');
%   sys = a.statespace(Modes = 20, Damping = 0.02);
%
% |statespace| returns A, B, C, D in modal coordinates [eta; deta], with
% signal names, the mode frequencies and mode shapes. Rigid-body and
% mechanism modes are double integrators. Joint and transmission dampers
% are added to the modal damping. Static correction (default)
% makes the DC gain exact with few modes. Use |Format = "ss"| for an |ss|
% object (Control System Toolbox), needed for |zpk|, |tf| and |feedback|;
% see <fem_analysis.html#4 Control design>.
%
%% Loads and outputs (by hand)
% Build loads and outputs in the full space and project them:
%
%   f  = zeros(a.nDof, 1);  f(d) = F;      % nodal forces
%   fq = T'*f;                             % reduced load
%   q  = K\fq;  u = T*q;                   % static solution
%   C  = zeros(1, a.nDof);  C(d(1)) = 1;  Cq = C*T;   % output
%
% Gravity is |T'*(Mfull*g)| with g = -9.81 in every z node dof. Support
% reactions are |Kfull*u - f| at the supported dofs. See
% <tut05_loadsStatics.html Tutorial 5>.
%
%% Plotting and animation
% |h = a.plot(u, scale)| draws all bodies deformed by the full vector u and
% returns one graphics handle per body. |a.animate(u, ...)| animates a mode
% in slow motion (vibration modes are far too fast for real time); the
% phase follows the wall clock, so the speed does not depend on rendering:
%
%   [V, f] = fem.modes(K, M, 6);
%   a.animate(T*V(:,4), Frequency = f(4), Title = 'Mode 4', Rate = 1, Duration = 3);
%
% Since u also holds rotation dofs (radians), scale with the node dofs
% only: |scale = s/max(abs(u(1:3*size(a.X, 1))))|.
