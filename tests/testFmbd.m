function testFmbd()
%TESTFMBD Verification of the fmbd package. Requires Gmsh (see fmbd.meshSTL).

addpath(fileparts(fileparts(mfilename('fullpath'))));
tmp = tempname;  mkdir(tmp);  cleanup = onCleanup(@() rmdir(tmp, 's'));
g = 9.81;  ode = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% double compound pendulum vs Lagrange equations (ode45)
m1 = 1.5;  m2 = 0.8;  L1 = 0.4;  lc1 = 0.2;  lc2 = 0.25;  I1 = 0.02;  I2 = 0.01;
s = fmbd.System(Gravity = [0 0 -g]);
s.add(fmbd.Body.pointMass(m1, [lc1 0 0], I1*eye(3), Name = 'link1'));
s.add(fmbd.Body.pointMass(m2, [L1+lc2 0 0], I2*eye(3), Name = 'link2'));
s.joint('ground', 'link1', 'revolute', Point = [0 0 0], Axis = [0 1 0], Name = 'j1');
s.joint('link1', 'link2', 'revolute', Point = [L1 0 0], Axis = [0 1 0], Name = 'j2');
s.output('ry1', 'link1', [], 5);  s.output('ry2', 'link2', [], 5);  s.output('q2', 'j2');
[t, y] = s.simulate([0 1], Step = 2.5e-4);
c12 = @(z) m2*L1*lc2*cos(z(1) - z(2));  s12 = @(z) m2*L1*lc2*sin(z(1) - z(2));
f = @(t, z) [z(3:4); [I1 + m1*lc1^2 + m2*L1^2, c12(z); c12(z), I2 + m2*lc2^2] \ ...
    -[s12(z)*z(4)^2 + (m1*lc1 + m2*L1)*g*sin(z(1)); -s12(z)*z(3)^2 + m2*lc2*g*sin(z(2))]];
[~, z] = ode45(f, t, [pi/2; pi/2; 0; 0], ode);
check('Double pendulum vs Lagrange (ode45)', pi/2 - unwrap([y.ry1 y.ry2]), z(:,1:2), 2e-5);
check('Joint angle == difference of link angles', wrapToPi(y.q2 - (y.ry2 - y.ry1)), 0, 1e-12);

% damped elastic belt at large angles == exact linear solution (parallel axes)
s = belt(true);  [t, y] = s.simulate([0 0.3], Step = 1e-4);
Kq = [50/9, -50/3; -50/3, 300 + 50];  Dq = [1e-3 + 0.02/9, -0.02/3; -0.02/3, 0.05 + 0.02];
Mq = diag([1e-4, 0.01 + 2*0.2^2]);  A = [zeros(2), eye(2); -Mq\Kq, -Mq\Dq];
X = cell2mat(arrayfun(@(tk) expm(A*tk)*[1.5; 0.5; 0; 0], t', 'UniformOutput', false))';
check('Belt at 1.5 rad == linear solution', [y.qm y.qa], X(:, 1:2), 2e-4);

% elastic transmission starts relaxed: the motor follows the link's Position
s = belt(false);  [~, y] = s.simulate([0 1e-3]);
check('Belt relaxed at start: motor = N x link Position', y.qm(1), 3*0.5, 1e-12);

% rigid transmission, motor torque, redundant support
N = 3;  Jm = 1e-4;  Js = 0.01 + 2*0.2^2;  tau = 0.02;
s = fmbd.System();
s.add(fmbd.Body.pointMass(2, [0.2 0 0], 0.01*eye(3), Name = 'link'));
s.add(fmbd.Body.pointMass(0.1, [0.5 0 0], Jm*eye(3), Name = 'rotor'));
s.joint('ground', 'link',  'revolute', Point = [0 0 0], Name = 'axis');
s.joint('ground', 'rotor', 'revolute', Point = [0.5 0 0], Name = 'motor');
s.transmission({'motor', 'axis'}, [1/N -1]);
s.fix('rotor', 1:3);                                     % redundant with the motor joint
s.input('tau', 'motor');  s.output('qm', 'motor');  s.output('qa', 'axis');
[t, y] = s.simulate([0 0.5], Step = 1e-3, Input = @(t, y) tau);
check('Rigid transmission + torque, (Jm + Js/N^2) qm'''' = tau', ...
      [y.qm, N*y.qa], repmat(0.5*tau/(Jm + Js/N^2)*t.^2, 1, 2), 1e-7);

% heavy top on a spherical joint, started by Position and Velocity
s = fmbd.System(Gravity = [0 0 -g]);
Jt = diag([2 2 4])*1e-3;
s.add(fmbd.Body.pointMass(1, [0 0 0.1], Jt, Name = 'top'));
s.joint('ground', 'top', 'spherical', Point = [0 0 0], Position = [0.3 0 0], ...
        Velocity = [0 0 60], Name = 'pivot');
[~, ~, x] = s.simulate([0 0.2], Step = 1e-4, Rho = 0.9);
E = s.energy(x);  R0 = x.R(:,:,1,1);  w0 = x.nu(4:6,1);  cb = R0'*x.r(:,1,1);
Jp = Jt + (cb'*cb*eye(3) - cb*cb');                    % body frame = world-aligned at definition
f = @(t, z) [reshape(reshape(z(1:9), 3, 3)*hat3(z(10:12)), 9, 1); ...
             Jp\(cross(cb, reshape(z(1:9), 3, 3)'*[0; 0; -g]) - cross(z(10:12), Jp*z(10:12)))];
[~, z] = ode45(f, [0 0.2], [R0(:); w0], ode);
check(sprintf('Heavy top (%.0f rad/s), energy', norm(w0)), E.total/E.total(1), 1, 1e-6);
check('Heavy top vs Euler equations (ode45)', x.r(:,1,end), reshape(z(end,1:9), 3, 3)*cb, 1e-6);

% free fall along a prismatic joint
s = fmbd.System(Gravity = [1 2 -g]);
s.add(fmbd.Body.pointMass(3, [0.1 0.2 0.3], 0.01*eye(3), Name = 'slider'));
s.joint('ground', 'slider', 'prismatic', Axis = [0 0 1], Name = 'z');
s.output('z', 'z');
[t, y] = s.simulate([0 0.3], Step = 1e-2);
check('Free fall on a prismatic joint', y.z, -0.5*g*t.^2, 1e-12);

% closed loop: Position on the crank of a four-bar == position analysis
O2 = [0 0];  A0 = [0 0.1];  B0 = [0.3 0.15];  O4 = [0.3 0];  th = 0.5;
s = fmbd.System();
pts = {O2, A0; A0, B0; B0, O4};  nm = {'crank', 'coupler', 'rocker'};
for i = 1:3
    c = ([pts{i,1}, 0] + [pts{i,2}, 0])/2;
    s.add(fmbd.Body.pointMass(0.2, c, 1e-3*eye(3), Name = nm{i}));
end
s.joint('ground',  'crank',   'revolute', Point = [O2 0], Position = th, Name = 'j1');
s.joint('crank',   'coupler', 'revolute', Point = [A0 0], Name = 'j2');
s.joint('coupler', 'rocker',  'revolute', Point = [B0 0], Name = 'j3');
s.joint('ground',  'rocker',  'revolute', Point = [O4 0], Name = 'j4');
s.output('q1', 'j1');  s.output('q4', 'j4');
[~, y] = s.simulate([0 1e-3]);
A1 = O2 + (A0 - O2)*[cos(th) sin(th); -sin(th) cos(th)];   % crank rotated by th (row vector)
b = norm(B0 - A0);  c = norm(B0 - O4);  d = norm(A1 - O4);
ang = atan2(A1(2) - O4(2), A1(1) - O4(1)) + [1 -1]*acos((c^2 + d^2 - b^2)/(2*c*d));
B1 = O4 + c*[cos(ang)', sin(ang)'];  [~, k] = min(vecnorm(B1 - B0, 2, 2));  B1 = B1(k,:);
q4 = atan2(B1(2) - O4(2), B1(1) - O4(1)) - atan2(B0(2) - O4(2), B0(1) - O4(1));
check('Four-bar: Position on the crank == position analysis', [y.q1(1), y.q4(1)], [th, q4], 1e-10);

% joints without a Position start at 0 (Position on the elbow of a chain)
s = fmbd.System();
s.add(fmbd.Body.pointMass(1, [0.2 0 0], 1e-3*eye(3), Name = 'l1'));
s.add(fmbd.Body.pointMass(1, [0.6 0 0], 1e-3*eye(3), Name = 'l2'));
s.joint('ground', 'l1', 'revolute', Point = [0 0 0], Name = 'shoulder');
s.joint('l1', 'l2', 'revolute', Point = [0.4 0 0], Position = pi/2, Name = 'elbow');
s.output('q1', 'shoulder');  s.output('q2', 'elbow');
[~, y] = s.simulate([0 1e-3]);
check('Position on the elbow only: shoulder stays at 0', [y.q1(1), y.q2(1)], [0, pi/2], 1e-12);

% operating-point inputs: an arm held at theta against its joint springs
% (equilibrium inputs, kept in x.u) == the same arm built at theta
th = 0.8;  kj = 200;
for built = [false true]
    s = fmbd.System();
    c2 = [0.2 0 0] + built*([cos(th) sin(th) 0]*0.1 - [0.1 0 0]) + [0.1 0 0];
    s.add(fmbd.Body.pointMass(0.1, [0.1 0 0], 1e-4*eye(3), Name = 'l1'));
    s.add(fmbd.Body.pointMass(0.1, c2, 1e-4*eye(3), Name = 'l2'));
    s.joint('ground', 'l1', 'revolute', Point = [0 0 0], Stiffness = kj, Name = 'shoulder');
    s.joint('l1', 'l2', 'revolute', Point = [0.2 0 0], Stiffness = kj, Name = 'elbow', ...
            Position = (~built)*th);
    s.input('tau1', 'shoulder');  s.input('tau2', 'elbow');
    if built, f2 = s.modes(); else, f1 = s.modes(s.equilibrium(Input = kj*[0; th])); end
end
check('Linearised with holding torques (x.u) == relaxed at that pose', f1, f2, 1e-9);

% linearisation: gravity stiffness of a hanging and an inverted pendulum
label = {'hanging pendulum: f = sqrt(m g l/J)', 'inverted pendulum: f = -sqrt(m g l/J)'};
for sg = [-1 1]                                          % centre of mass below/above the hinge
    s = fmbd.System(Gravity = [0 0 -g]);
    s.add(fmbd.Body.pointMass(2, [0 0 sg*0.3], 0.01*eye(3), Name = 'p'));
    s.joint('ground', 'p', 'revolute', Point = [0 0 0], Axis = [0 1 0], Name = 'h');
    check(['Linearised ' label{(sg + 3)/2}], s.modes(), -sg*sqrt(2*g*0.3/(0.01 + 2*0.3^2))/(2*pi), 1e-10);
end

% equilibrium of a pendulum on a torsion spring, k th = m g l cos(th), and
% the linearised frequency there
s = fmbd.System(Gravity = [0 0 -g]);
s.add(fmbd.Body.pointMass(2, [0.3 0 0], 0.01*eye(3), Name = 'p'));
s.joint('ground', 'p', 'revolute', Point = [0 0 0], Axis = [0 1 0], Stiffness = 10, Name = 'h');
s.output('th', 'h');
x = s.equilibrium();  [~, y] = s.simulate([0 1e-3], Initial = x);
th = fzero(@(th) 10*th - 2*g*0.3*cos(th), 0.5);
check('Equilibrium: k th = m g l cos(th)', y.th(1), th, 1e-10);
check('Linearised at the equilibrium (preloaded spring)', s.modes(x), ...
      sqrt((10 + 2*g*0.3*sin(th))/(0.01 + 2*0.3^2))/(2*pi), 1e-8);

% equilibrium: converges on a gravity-loaded four-bar 1 km from the origin
% (relative constraint tolerance), and stops early when there is none
P = {[0 0 0], [0.05 0.1 0], [0.3 0.15 0], [0.25 0 0]};  O = [1e3 0 0];
s = fmbd.System(Gravity = [0 -g 0]);
for i = 1:3
    s.add(fmbd.Body.pointMass(0.3, O + (P{i} + P{i+1})/2, 1e-3*eye(3), Name = sprintf('l%d', i)));
end
s.joint('ground', 'l1', 'revolute', Point = O + P{1}, Stiffness = 5);
s.joint('l1', 'l2', 'revolute', Point = O + P{2});
s.joint('l2', 'l3', 'revolute', Point = O + P{3});
s.joint('ground', 'l3', 'revolute', Point = O + P{4});
[~, info] = s.equilibrium();
check('Equilibrium far from the origin converges (<= 10 iterations)', ...
      [info.converged, info.iterations <= 10], [1 1], 0);
s = fmbd.System(Gravity = [0 0 -g]);                    % free horizontal pendulum: no equilibrium
s.add(fmbd.Body.pointMass(1, [0.3 0 0], 1e-3*eye(3), Name = 'p'));
s.joint('ground', 'p', 'revolute', Point = [0 0 0], Axis = [0 1 0]);
w = warning('off', 'fmbd:equilibrium');  [~, info] = s.equilibrium();  warning(w);
check('No equilibrium: stops early and reports it', [info.converged, info.iterations <= 10], [0 1], 0);

% equilibrium with stiff elements and a free mechanism: an elastic plate
% bolted at two couplings to a link on a free joint, driven by a stiff belt,
% under gravity (constraints converge; the free angles stay at zero)
plate = fullfile(tmp, 'plate.stl');  fmbd.boxSTL(plate, [0.2 -0.05 -0.003], [0.3 0.05 0.003]);
s = fmbd.System(Gravity = [0 0 -g]);
s.add(fmbd.Body.pointMass(1, [0.1 0 0], 1e-3*eye(3), Name = 'link'));
s.add(fmbd.Body.pointMass(0.1, [0 0.1 0], 0.1*eye(3), Name = 'rotor'));
s.add(fmbd.Body.elastic(plate, fmbd.material(70e9, 0.33, 2700), Size = 0.01, Order = 1, ...
                        Modes = 10, Name = 'plate'));
s.joint('ground', 'link', 'revolute', Point = [0 0 0], Stiffness = [1e5 1e5 0], ...
        Free = [0 0 0 1 1 1], Name = 'shoulder');
s.joint('ground', 'rotor', 'revolute', Point = [0 0.1 0], Name = 'motor');
s.transmission({'motor', 'link'}, [1 -1], Dirs = {1, 6}, Stiffness = 1.7e5);
for k = 1:2
    h = s.couple('plate', @(X) X(:,1) < 0.21 & (2*k - 3)*X(:,2) > 0.03, Name = sprintf('bolt%d', k));
    s.joint('link', h, 'fixed');
end
s.output('m', 'motor');  s.output('rz', 'link', [], 6);
[x, info] = s.equilibrium();  [~, y] = s.simulate([0 1e-4], Initial = x, Step = 1e-4);
check('Equilibrium: stiff belt, bolted plate, free mechanism (converged, angles 0)', ...
      [~info.converged, info.constraint, y.m(1), y.rz(1)], 0, 1e-9);

% state space of the belt drive == the exact frequency response, and the DC gain
% with one mode plus static correction
s = belt(false);  s.input('tau', 'motor');  s.output('wa', 'axis', Type = 'velocity');
S1 = s.statespace();  S3 = s.statespace(Modes = 1);
Hs = @(S, w) S.C*((1i*w*eye(size(S.A)) - S.A)\S.B) + S.D;
Kq = [50/9, -50/3; -50/3, 300 + 50];  Dq = [1e-3 + 0.02/9, -0.02/3; -0.02/3, 0.05 + 0.02];
Mq = diag([1e-4, 0.01 + 2*0.2^2]);
Hx = @(w) (-w^2*Mq + 1i*w*Dq + Kq)\[1; 0];                          % [motor; axis] per unit torque
w = 2*pi*[0.1 3 10 60 300];  err = 0;
for k = 1:numel(w)
    Hq = Hx(w(k));  Hr = [Hq(1); Hq(2); 1i*w(k)*Hq(2)];             % qm, qa, wa
    err = max(err, norm(Hs(S1, w(k)) - Hr)/norm(Hr));
end
check('State space: belt drive FRF == exact', err, 0, 1e-9);
check('State space: 1 mode + static correction keeps the DC gain', ...
      S3.C(1,:)*(S3.A\S3.B), S1.C(1,:)*(S1.A\S1.B), 1e-10);

% elastic bodies: cantilever with a tip mass on its second coupling
Lb = 0.5;  hb = 0.02;  beam = fullfile(tmp, 'beam.stl');
fmbd.boxSTL(beam, [0 0 0], [Lb hb hb]);
steel = fmbd.material(210e9, 0.3, 7800);
[~, y] = elasticBeamTip(beam, steel, Lb, hb, 10).simulate([0 0.02], Step = 4e-5, ...
                                                       Input = @(t, y) 1, Rho = 0.95);
[~, yr] = elasticBeamTip(beam, steel, Lb, hb, 40).simulate([0 0.02], Step = 4e-5, ...
                                                        Input = @(t, y) 1, Rho = 0.95);
check('Elastic: tip mass and node set, 10 modes == 40 modes', [y.uz, y.mid]./max(abs([yr.uz yr.mid])), ...
      [yr.uz yr.mid]./max(abs([yr.uz yr.mid])), 1e-3);

% equilibrium under gravity and the first mode of a Tet10 cantilever
% (clamped on nodes) == Euler-Bernoulli
s = fmbd.System(Gravity = [0 0 -g]);
s.add(fmbd.Body.elastic(beam, steel, Size = hb/2, Name = 'beam', Modes = 4));
s.fix('beam', @(X) X(:,1) < 1e-6);
s.output('tipz', 'beam', @(X) X(:,1) > Lb - 1e-6, 3);
xe = s.equilibrium();  [~, ye] = s.simulate([0 1e-4], Initial = xe, Step = 1e-4);
EI = 210e9*hb^4/12;  q = 7800*hb^2*g;
check('Elastic: gravity sag at equilibrium == q L^4/(8 EI)', ye.tipz(1)/(-q*Lb^4/(8*EI)), 1, 1e-2);
check('Elastic: first mode at equilibrium == Euler-Bernoulli', s.modes(xe, 1), ...
      1.875104^2/(2*pi)*sqrt(EI/(7800*hb^2*Lb^4)), 1e-2);

% elastic body tied to a rigid body (a coupling + fixed joint) == tied to a
% very stiff elastic body
b3 = {fullfile(tmp, 'l3.stl'), fullfile(tmp, 'r3.stl')};
fmbd.boxSTL(b3{1}, [0 0 0], [Lb/2 hb hb]);  fmbd.boxSTL(b3{2}, [Lb/2 0 0], [Lb hb hb]);
f = cell(1, 2);
for rigid = [true false]
    s = fmbd.System();
    s.add(fmbd.Body.elastic(b3{1}, steel, Size = hb/2, Name = 'left', Modes = 20));
    if rigid
        s.add(fmbd.Body(b3{2}, steel, Size = hb/2, Name = 'right', Rigid = true));
    else
        stiff = fmbd.material(210e15, 0.3, 7800);
        s.add(fmbd.Body.elastic(b3{2}, stiff, Size = hb/2, Name = 'right', Modes = 4));
    end
    s.fix('left', @(X) X(:,1) < 1e-6);  s.tie('right', 'left');
    f{rigid + 1} = s.modes([], 4);
end
check('Tie elastic-rigid == tie to a very stiff elastic body', f{2}, f{1}, 5e-4);

% acceleration outputs: mass on a spring under gravity and a step force
s = fmbd.System(Gravity = [0 0 -g]);
s.add(fmbd.Body.pointMass(2, [0 0 0], 1e-3*eye(3), Name = 'm'));
s.joint('ground', 'm', 'prismatic', Axis = [0 0 1], Stiffness = 800, Name = 'z');
s.input('F', 'z');  s.output('a', 'z', Type = 'acceleration');
s.output('az', 'm', [], 3, Type = 'acceleration');
[t, y] = s.simulate([0 0.2], Step = 1e-4, Input = @(t, y) 3);
check('Acceleration outputs (joint, point) vs analytic', [y.a y.az], repmat((1.5 - g)*cos(20*t), 1, 2), 2e-5);
S = s.statespace();
check('State space: acceleration feedthrough D = 1/m', S.D, [0.5; 0.5], 1e-12);

% ties and node supports: non-matching halves, clamped on nodes, vs the
% monolithic beam; a stiff elastic interface vs the rigid tie
half = {fullfile(tmp, 'left.stl'), fullfile(tmp, 'right.stl'), fullfile(tmp, 'full.stl')};
fmbd.boxSTL(half{1}, [0 0 0], [Lb/2 hb hb]);  fmbd.boxSTL(half{2}, [Lb/2 0 0], [Lb hb hb]);
fmbd.boxSTL(half{3}, [0 0 0], [Lb hb hb]);
clamp = @(X) X(:,1) < 1e-6;
s = fmbd.System();
s.add(fmbd.Body.elastic(half{3}, steel, Size = hb/2, Name = 'beam', Modes = 20));
s.fix('beam', clamp);  fmono = s.modes([], 4);
for kind = ["rigid", "elastic"]
    s = fmbd.System();
    s.add(fmbd.Body.elastic(half{1}, steel, Size = hb/2, Name = 'left', Modes = 10));
    s.add(fmbd.Body.elastic(half{2}, steel, Size = hb/3, Name = 'right', Modes = 10));
    s.fix('left', clamp);
    if kind == "rigid", ns = s.tie('left', 'right'); else, s.tie('left', 'right', Stiffness = 1e16); end
    check(sprintf('Tie (%s, %d slave nodes): halves == monolithic beam', kind, ns), ...
          s.modes([], 4), fmono, 1e-3);
end

% static correction: a node-set load away from couplings keeps the exact
% static gain even without vibration modes
s = fmbd.System();
s.add(fmbd.Body.elastic(half{3}, steel, Size = hb/2, Name = 'beam', Modes = 0));
s.fix('beam', clamp);
tipSel = @(X) X(:,1) > Lb - 1e-6;
s.input('F', 'beam', tipSel, 3);  s.output('w', 'beam', tipSel, 3);
S = s.statespace();  dc = -S.C*(S.A\S.B);
s.Bodies{1}.Modes = 30;  s.update();             % many modes, same load shape
S30 = s.statespace();
check('Static correction: node-set load, DC gain with 0 modes == 30 modes', dc, -S30.C*(S30.A\S30.B), 1e-9);

% elastic bodies: stiff flexible pendulum == rigid pendulum (3 rad swing)
stiff = fmbd.material(210e12, 0.3, 7800);
[~, yr] = pendulum(beam, stiff, Lb, true, g).simulate([0 0.6]);
[~, ye] = pendulum(beam, stiff, Lb, false, g).simulate([0 0.6]);
check(sprintf('Elastic: stiff pendulum (%.1f rad) == rigid', max(abs(yr.th))), ye.th, yr.th, 1e-4);

% elastic bodies: energy conservation; kinetic energy from the inertia
% invariants == nodal 1/2 xdot' M xdot
s = pendulum(beam, fmbd.material(2e9, 0.3, 7800), Lb, false, g);
[~, ~, x] = s.simulate([0 0.6], Step = 5e-4, Rho = 1);
E = s.energy(x);  k = size(x.nu, 2);  [Phi, X0] = s.modeShapes('beam');
[~, Mb] = s.Bodies{1}.matrices();
U = X0' + reshape(Phi*x.e(:,k), 3, []);  w = x.nu(4:6,k);
xd = x.nu(1:3,k) + x.R(:,:,1,k)*(cross(repmat(w, 1, size(U, 2)), U) + ...
     reshape(Phi*x.nu(7:end,k), 3, []));
check(sprintf('Elastic: soft pendulum energy (strain %.0e J)', max(E.strain)), ...
      (E.total - E.total(1))/max(E.kinetic), 0, 1e-5);
check('Elastic: kinetic energy, invariants == nodal', E.kinetic(k), 0.5*full(xd(:)'*Mb*xd(:)), 1e-10);
end

function s = belt(posed)
% Link on a torsion spring, driven by a rotor through a damped elastic belt.
s = fmbd.System();
s.add(fmbd.Body.pointMass(2, [0.2 0 0], 0.01*eye(3), Name = 'link'));
s.add(fmbd.Body.pointMass(0.1, [0.5 0 0], 1e-4*eye(3), Name = 'rotor'));
pm = [];  if posed, pm = 1.5; end
s.joint('ground', 'link',  'revolute', Point = [0 0 0], Stiffness = 300, Damping = 0.05, ...
        Position = 0.5, Name = 'axis');
s.joint('ground', 'rotor', 'revolute', Point = [0.5 0 0], Damping = 1e-3, Position = pm, Name = 'motor');
s.transmission({'motor', 'axis'}, [1/3 -1], Stiffness = 50, Damping = 0.02);
s.output('qm', 'motor');  s.output('qa', 'axis');
end

function s = elasticBeamTip(file, mat, L, h, modes)
% Clamped elastic beam with a point mass on its tip coupling.
s = elasticBeam(file, mat, L, [0 2e-5], modes);
s.add(fmbd.Body.pointMass(0.3, [L + 0.03, h/2, h/2], 1e-4*eye(3), Name = 'mass'));
s.joint('ground', 'root', 'fixed');  s.joint('tip', 'mass', 'fixed', Point = [L h/2 h/2]);
s.input('Fz', 'mass', [], 3);  s.output('uz', 'mass', [], 3);
s.output('mid', 'beam', @(X) abs(X(:,1) - L/2) < 0.006, 3);
end

function s = elasticBeam(file, mat, L, damping, modes)
% Elastic Tet4 beam along x with couplings at both ends.
s = fmbd.System();
s.add(fmbd.Body.elastic(file, mat, Size = 0.01, Order = 1, Name = 'beam', Damping = damping, Modes = modes));
s.couple('beam', @(X) X(:,1) < 1e-9, Name = 'root');
s.couple('beam', @(X) X(:,1) > L - 1e-9, Name = 'tip');
end

function s = pendulum(file, mat, L, rigid, g)
% Beam pendulum hinged about y at its root, rigid or elastic.
if rigid
    s = fmbd.System();
    s.add(fmbd.Body(file, mat, Size = 0.01, Order = 1, Name = 'beam', Rigid = true));
    s.joint('ground', 'beam', 'revolute', Point = [0 0.01 0.01], Axis = [0 1 0], Name = 'hinge');
else
    s = elasticBeam(file, mat, L, [0 0], 6);
    s.joint('ground', 'root', 'revolute', Point = [0 0.01 0.01], Axis = [0 1 0], Name = 'hinge');
end
s.Gravity = [0 0 -g];
s.output('th', 'hinge');
end

function S = hat3(v), S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]; end

function check(name, val, ref, tol)
err = max(abs(val(:) - ref(:))./max(abs(ref(:)), 1));
if isequal(ref, 0), err = max(abs(val(:))); end
status = "PASS";  if err > tol, status = "FAIL"; end
fprintf('%-4s  %-55s err = %.2e\n', status, name, err);
end
