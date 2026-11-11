%% Tutorial 1 - Your first model: a two-link arm with an elastic forearm
% This tutorial goes through the basic workflow:
%
% # geometry (STL) -> |fembd.Body|: rigid, or elastic (Gmsh mesh + material)
% # bodies -> |fembd.Assembly|: joints, couplings, inputs and outputs
% # |a.nonlinear()| -> |simulate|: large motion with a controller in the loop
% # |a.system()|, |a.statespace()|: linear analysis about a pose
%
% The arm moves in the horizontal plane (a SCARA arm): a rigid upper arm on
% a shoulder joint, and an elastic aluminium forearm on an elbow joint, both
% driven by a motor torque under PD control. Gravity bends the forearm.
% Run it section by section (Ctrl+Enter), or render it with
% |publish('tut01_firstModel')|.

addpath(fileparts(fileparts(mfilename('fullpath'))));     % toolbox root
tmp = tempname;  mkdir(tmp);

%% 1. Geometry and bodies
% Two boxes as STL files (in m). The upper arm is rigid: its mesh only gives
% its mass, inertia and surface. The forearm is elastic: Gmsh meshes it with
% quadratic Tet10 elements.
fembd.boxSTL(fullfile(tmp, 'upper.stl'),   [0 -0.02 -0.01],     [0.3 0.02 0.01]);
fembd.boxSTL(fullfile(tmp, 'forearm.stl'), [0.3 -0.015 -0.005], [0.6 0.015 0.005]);
steel = fembd.material(210e9, 0.30, 7850);
alu   = fembd.material(70e9,  0.33, 2700);

upper   = fembd.Body(fullfile(tmp, 'upper.stl'),   steel, Size = 0.02, Rigid = true, Name = 'upper');
forearm = fembd.Body(fullfile(tmp, 'forearm.stl'), alu,   Size = 0.006, Name = 'forearm');
fprintf('upper arm: %.2f kg (rigid); forearm: %.3f kg, %d nodes\n', ...
        upper.massProperties(), forearm.massProperties(), forearm.nNode);

%% 2. Assembly: joints, couplings, inputs and outputs
% Joints connect reference points: a rigid body, or a patch of an elastic
% body made rigid by a coupling. The forearm gets couplings at its root
% (the elbow) and at its tip (the payload point).
a = fembd.Assembly();
a.add(upper);  a.add(forearm);
a.couple('forearm', @(X) X(:,1) < 0.3 + 1e-6, Name = 'root', Point = [0.3 0 0]);
a.couple('forearm', @(X) X(:,1) > 0.6 - 1e-6, Name = 'tip');
a.joint('ground', 'upper', 'revolute', Point = [0 0 0],   Name = 'shoulder');   % about z
a.joint('upper',  'root',  'revolute', Point = [0.3 0 0], Name = 'elbow');

% motor torques on the joints, and sensors: joint angles and rates, and the
% vertical tip motion (the forearm's bending)
a.input('tau1', 'shoulder');  a.input('tau2', 'elbow');
a.output('th1', 'shoulder');  a.output('w1', 'shoulder', Type = 'velocity');
a.output('th2', 'elbow');     a.output('w2', 'elbow',    Type = 'velocity');
a.output('tipz', 'tip', [], 3);

%% 3. Nonlinear simulation with a controller in the loop
% |a.nonlinear| builds the nonlinear model: exact joint kinematics for large
% angles, and the forearm's deformation by a few Craig-Bampton modes in a
% frame that moves with it. |simulate| integrates in time; the controller
% |Input = @(t, y)| gets the outputs of the previous step.
m = a.nonlinear(Gravity = [0 0 -9.81], Modes = 6);

ref = [1.0; -1.5];                                     % target angles (rad)
kp  = [16; 1.5];  kd = [1.8; 0.15];                    % PD gains (about 2 Hz bandwidth)
pd  = @(t, y) kp.*(ref - [y.th1; y.th2]) - kd.*[y.w1; y.w2];
[t, y, x] = m.simulate([0 1.5], Step = 5e-4, Input = pd);

figure('Color', 'w');
subplot(2, 1, 1);  plot(t, [y.th1 y.th2], 'LineWidth', 1.5);  grid on
ylabel('angle (rad)');  legend('shoulder', 'elbow', 'Location', 'east');
subplot(2, 1, 2);  plot(t, 1e6*y.tipz, 'LineWidth', 1.5);  grid on
xlabel('time (s)');  ylabel('tip z (\mum)');  title('forearm bending: gravity sag and vibration');
fprintf('final angles %.3f, %.3f rad (targets %.1f, %.1f)\n', y.th1(end), y.th2(end), ref);

%% 4. Energy and animation
% |energy| splits the energy into kinetic, gravity, strain and spring
% parts. |animate| plays the motion in real time (|Rate|); |Scale|
% magnifies the deformation.
E = m.energy(x);
fprintf('largest strain energy in the forearm: %.2e J\n', max(E.strain));
figure('Color', 'w');
m.animate(t, x, Rate = 1, Scale = 200, Title = 'SCARA arm');

%% 5. Linear analysis about a pose
% The same assembly gives the linear finite-element model about its current
% pose: all constraints are eliminated with u = T*q, and K and M are the
% stiffness and mass in the independent dofs q. The two joint angles are
% free, so the first two modes are rigid-body (mechanism) modes at 0 Hz.
[K, M, T] = a.system();
[V, f] = fembd.modes(K, M, 6);
disp(table((1:6)', f, 'VariableNames', {'mode', 'f_Hz'}));

% the state-space model from the motor torques to the sensors (for control
% design, see Tutorial 9 and the manual)
sys = a.statespace(Modes = 10, Damping = 0.01);
fprintf('state-space model: %d states, inputs %s, outputs %s\n', size(sys.A, 1), ...
        strjoin(sys.InputName, ', '), strjoin(sys.OutputName, ', '));
