%% Modal Analysis, Dynamics and Reduction
%
%% Natural frequencies
%
%   [V, f] = fem.modes(K, M, n)
%
% The lowest |n| frequencies (Hz) and mass-normalised mode shapes
% (V'*M*V = I), by shift-invert Lanczos. Free and partly free systems
% (singular K) are handled: rigid-body and mechanism modes appear with
% f = 0. Small systems are solved densely.
%
%% Dynamics
% The equations of motion in q are
%
% $$M\ddot q + D\dot q + K q = T^\top f(t), \qquad u = T q .$$
%
% fem returns K and M, and D for body, joint and transmission dampers
% (|[K, M, T, D] = a.system()|). Material damping of a body is Rayleigh
% damping, D = alpha*M + beta*K, set with |fem.Body(..., Damping = [alpha beta])|;
% |beta = 2*zeta/omega| gives the ratio zeta at omega. <tut08_dynamics.html Tutorial 8> shows frequency
% responses (direct and modal), Newmark time integration, and state-space
% forms. For explicit time integration use |a.system(Lumped = true)|.
%
%% State-space models
% With inputs and outputs defined (|a.input|, |a.output|), one call gives a
% reduced modal state-space model, including static correction, damping
% and acceleration feedthrough:
%
%   sys = a.statespace(Modes = 20, Damping = 0.02);
%
% See <fem_system.html Dofs, Supports, Loads and the System>.
%
%% Control design
% With the Control System Toolbox, ask for an |ss| object: |zpk|, |tf|,
% |bode| and |feedback| do not accept the default struct (the error is
% "Invalid syntax"). Channels are selected by signal name:
%
%   sys = a.statespace(Modes = 20, Format = "ss");
%   s   = tf('s');
%   G11 = zpk(sys('enc1', 'tau1'));              % one channel
%   G11r = tf(balred(sys('enc1', 'tau1'), 6));   % low-order tf in s
%
% A decentralised PID with a second-order low-pass per loop, closed around
% the full 2 x 2 plant so the coupling between the loops is kept:
%
%   P = sys({'enc1', 'enc2'}, {'tau1', 'tau2'});
%   C = tf(zeros(2));
%   for k = 1:2
%       LP     = wlp(k)^2/(s^2 + 2*zlp(k)*wlp(k)*s + wlp(k)^2);
%       C(k,k) = (KP(k) + KI(k)/s + KD(k)*s)*LP;
%   end
%   CL = feedback(P*C, eye(2));                  % closed loop
%   margin(C(1,1)*P(1,1))                        % loop-at-a-time margins
%
% Notes:
%
% * Keep the full model as |ss|: a |tf| of 40+ states has badly
%   conditioned coefficients. Use |zpk| per channel, or reduce with
%   |balred| before converting to |tf|.
% * A pure derivative |KD*s| is improper; roll it off with a low-pass (as
%   above) or use |pid(KP, KI, KD, Tf)|.
% * Torque to its own joint or motor coordinate (|tau1 -> enc1|) is
%   collocated and robust. A coordinate behind an elastic transmission
%   (|tau1 -> th1|, the link side of a belt) is not: the transmission mode
%   adds a 180 degree phase drop and limits the bandwidth.
% * With |Damping = 0| only the physical dampers act; a small modal ratio
%   (0.01-0.02) gives more realistic margins.
%
%% Craig-Bampton reduction
%
%   [Kc, Mc, V, f] = fem.craigBampton(K, M, b, n)
%
% Keeps the dofs |b| (indices into K, e.g. from |a.qdof|) as physical
% coordinates and adds |n| fixed-interface modes:
%
% $$q \approx V \left[\begin{array}{c} u_b \\ \eta \end{array}\right], \quad
%   V = \left[\begin{array}{cc} I & 0 \\ \Psi & \Phi \end{array}\right], \quad
%   \Psi = -K_{ii}^{-1} K_{ib}.$$
%
% |n = 0| gives Guyan reduction, which is exact for static loads at the
% boundary. The full field follows from u = T*V*x. Good boundary dofs are
% interfaces, actuators and sensors: e.g. the 6 dofs of a coupling or
% rigid body, or joint coordinates. See <tut09_reduction.html Tutorial 9>.
%
%   ib = a.qdof(a.dof('flange'));
%   [Kc, Mc, Vc] = fem.craigBampton(K, M, ib, 20);
%
% The 6 x 6 stiffness at a rigid body follows from Guyan reduction onto its
% dofs: |Ks = fem.craigBampton(K, M, a.qdof(a.dof('shuttle')), 0)|.
