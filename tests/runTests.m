function runTests()
%RUNTESTS Verification of the fembd package. Requires Gmsh (see fembd.meshSTL).

addpath(fileparts(fileparts(mfilename('fullpath'))));
tmp = tempname;  mkdir(tmp);  cleanup = onCleanup(@() rmdir(tmp, 's'));
steel = fembd.material(210e9, 0.3, 7800);
Lx = 1;  h = 0.05;

% single body: mass, rigid-body invariance, node ordering
fembd.boxSTL(fullfile(tmp, 'beam.stl'), [0 0 0], [Lx h h]);
for order = [1 2]
    b = fembd.Body(fullfile(tmp, 'beam.stl'), steel, Size = h/2, Order = order);
    [K, M] = b.matrices();
    e = repmat([1; 0; 0], b.nNode, 1);
    check(sprintf('Tet%d mass', 4 + 6*(order == 2)), e'*M*e, steel.rho*Lx*h*h, 1e-6);   % STL is single precision
    check(sprintf('Tet%d rigid translation', 4 + 6*(order == 2)), norm(K*e)/norm(K, 1), 0, 1e-10);
    r = reshape([zeros(1, b.nNode); -b.X(:,3)'; b.X(:,2)'], [], 1);   % rotation about x
    check(sprintf('Tet%d rigid rotation', 4 + 6*(order == 2)), norm(K*r)/norm(K, 1)/norm(r), 0, 1e-10);
end
mid = b.X(b.T(:, 5:10), :);
cor = (b.X(b.T(:, [1 2 3 4 4 4]), :) + b.X(b.T(:, [2 3 1 1 3 2]), :))/2;
check('Tet10 edge-node ordering', max(vecnorm(mid - cor, 2, 2)), 0, 1e-12);

% cantilever: first bending frequency vs Euler-Bernoulli
a = fembd.Assembly();  a.add(b);
a.fix(1, @(X) X(:,1) < 1e-9);
[Kr, Mr] = a.system();
[~, f] = fembd.modes(Kr, Mr, 4);
fEB = 1.875104^2/(2*pi)*sqrt(210e9*h^4/12/(7800*h^2*Lx^4));
check('Cantilever f1 vs Euler-Bernoulli', f(1), fEB, 1e-2);
fRef = f;

% lumped (HRZ) mass
[~, Ml] = b.matrices(Lumped = true);
check('Lumped mass: diagonal, positive', nnz(Ml - diag(diag(Ml))) + any(diag(Ml) <= 0), 0, 0);
check('Lumped mass: total', e'*Ml*e, steel.rho*Lx*h*h, 1e-6);
[Kr, Mr] = a.system(Lumped = true);
[~, f] = fembd.modes(Kr, Mr, 4);
check('Lumped cantilever vs consistent', f, fRef, 2e-2);

% free body: six rigid-body modes
a = fembd.Assembly();  a.add(b);
[Kr, Mr] = a.system();
[~, f] = fembd.modes(Kr, Mr, 7);
check('Free body: 6 rigid modes', max(f(1:6))/f(7), 0, 1e-4);

% Rayleigh body damping: modal ratios alpha/(2 w) + beta w/2
aR = 2;  bR = 2e-5;
a = fembd.Assembly();
a.add(fembd.Body(fullfile(tmp, 'beam.stl'), steel, Size = h/2, Order = 1, Damping = [aR bR]));
a.fix(1, @(X) X(:,1) < 1e-9);
a.input('F', 1, @(X) X(:,1) > Lx - 1e-9, 3);  a.output('uz', 1, @(X) X(:,1) > Lx - 1e-9, 3);
[Kr, Mr, ~, Dr] = a.system();
check('Body damping: D = alpha M + beta K', norm(Dr - aR*Mr - bR*Kr, 1)/norm(Dr, 1), 0, 1e-12);
sys = a.statespace(Modes = 4, Damping = 0, StaticCorrection = false);
p = eig(sys.A);  [~, i] = sort(abs(p));  p = p(i(1:2:end));
w = 2*pi*sys.f;
check('Body damping: modal ratios alpha/(2w) + beta w/2', -real(p)./abs(p), aR./(2*w) + bR*w/2, 1e-8);
ok = false;  try, fembd.Body(struct('X', b.X, 'T', b.T), steel, Rigid = true, Damping = [1 0]); catch, ok = true; end
check('Body damping on a rigid body is rejected', ok, true, 0);

% two non-conforming halves, tied
fembd.boxSTL(fullfile(tmp, 'left.stl'),  [0 0 0],    [Lx/2 h h]);
fembd.boxSTL(fullfile(tmp, 'right.stl'), [Lx/2 0 0], [Lx h h]);
L = fembd.Body(fullfile(tmp, 'left.stl'),  steel, Size = h/2);
R = fembd.Body(fullfile(tmp, 'right.stl'), steel, Size = h/3);
a = fembd.Assembly();  a.add(L);  a.add(R);
ns = a.tie('left', 'right');
a.fix('left', @(X) X(:,1) < 1e-9);
[Kr, Mr, T] = a.system();
[V, f] = fembd.modes(Kr, Mr, 4);
check(sprintf('Tied halves (%d slave nodes) vs monolithic', ns), f, fRef, 5e-3);
U = reshape(T*V(:,1), 3, [])';  lk = a.Links(1);
Um = reshape(U(lk.master, :), [size(lk.master), 3]);
gap = U(lk.slave, :) - squeeze(sum(lk.N.*Um, 2));
check('Tie: interface gap / max displacement', max(abs(gap(:)))/max(abs(U(:))), 0, 1e-12);

% Craig-Bampton: tip face as boundary, 20 fixed-interface modes
ib = a.qdof(a.dof('right', @(X) X(:,1) > Lx - 1e-9));
[Kc, Mc] = fembd.craigBampton(Kr, Mr, ib, 20);
[~, fc] = fembd.modes(Kc, Mc, 4);
check(sprintf('Craig-Bampton (%d+20 of %d dofs) vs full', numel(ib), size(Kr, 1)), fc, f, 1e-3);
check('Craig-Bampton: K boundary/modal decoupled', ...
      norm(Kc(1:numel(ib), numel(ib)+1:end))/norm(Kc), 0, 1e-8);

% elastic interface with very high stiffness approaches the tie
a.Links = a.Links([]);
a.tie('left', 'right', Stiffness = 1e18*[1 1]);
check('Spring interface area', sum(a.Links(1).A), h^2, 1e-6);
[Kr, Mr] = a.system();
[~, f] = fembd.modes(Kr, Mr, 4);
check('Stiff spring interface vs tie', f, fRef, 5e-3);

% rigid bodies: exact 6x6 mass / inertia of a free rigid box
Rb = fembd.Body(fullfile(tmp, 'beam.stl'), struct('rho', 7800), Size = h, Order = 1, Rigid = true);
a = fembd.Assembly();  a.add(Rb);
[Kr, Mr] = a.system();
m = 7800*Lx*h*h;
Mref = diag(m*[1 1 1 (2*h^2)/12 (Lx^2 + h^2)/12 (Lx^2 + h^2)/12]);
check('Rigid box: K = 0, M = rigid mass/inertia', ...
      [norm(Kr, 1), norm(full(Mr) - Mref)/m], [0 0], 1e-6);

% beam with a rigid right half (slave given first: swapped to master)
% vs a nearly rigid elastic right half
Rr = fembd.Body(fullfile(tmp, 'right.stl'), steel, Size = h/2, Rigid = true);
Rs = fembd.Body(fullfile(tmp, 'right.stl'), fembd.material(210e9*1e5, 0.3, 7800), Size = h/2);
right = {Rr, Rs};
for k = 1:2
    a = fembd.Assembly();  a.add(L);  a.add(right{k});
    a.tie('left', 2);
    a.fix('left', @(X) X(:,1) < 1e-9);
    [Kr, Mr] = a.system();
    [~, fk(:,k)] = fembd.modes(Kr, Mr, 4);   %#ok<AGROW>
end
check(sprintf('Rigid half (%d dofs) vs E x 1e5 half', size(Kr, 1)), fk(:,1), fk(:,2), 1e-3);

% rigid body fixed via its reference dofs == clamping the tied face
a = fembd.Assembly();  a.add(L);  a.add(Rr);  a.tie(2, 1);
a.fix(1, @(X) X(:,1) < 1e-9);  a.fix(2, [], 1:6);
[Kr, Mr] = a.system();  [~, f1] = fembd.modes(Kr, Mr, 4);
a = fembd.Assembly();  a.add(L);
a.fix(1, @(X) X(:,1) < 1e-9 | X(:,1) > Lx/2 - 1e-9);
[Kr, Mr] = a.system();  [~, f2] = fembd.modes(Kr, Mr, 4);
check('Fixed rigid body == clamped interface', f1, f2, 1e-10);
fRig = fk(:,1);

% mesh-free rigid bodies: mass properties from the STL surface
Sb = fembd.Body.rigidSurface(fullfile(tmp, 'beam.stl'), 7800);
check('Rigid surface: mass and inertia of a box', [Sb.MassProps.m, diag(Sb.MassProps.J)'], ...
      [m, diag(Mref(4:6, 4:6))'], 1e-6);
fembd.cylinderSTL(fullfile(tmp, 'cyl.stl'), [0 0], 0.1, 0, 0.5, 32);
Sc = fembd.Body.rigidSurface(fullfile(tmp, 'cyl.stl'), 1);
check('Rigid surface: volume of a 32-gon prism', Sc.MassProps.m, 32/2*0.1^2*sin(2*pi/32)*0.5, 1e-6);

% the rigid right half as surface only, and as a point mass + joint
mR = 7800*(Lx/2)*h*h;  JR = mR/12*diag([2*h^2, (Lx/2)^2 + h^2, (Lx/2)^2 + h^2]);
a = fembd.Assembly();  a.add(L);  a.add(fembd.Body.rigidSurface(fullfile(tmp, 'right.stl'), 7800));
a.tie('left', 'right');  a.fix('left', @(X) X(:,1) < 1e-9);
[Kr, Mr] = a.system();  [~, f] = fembd.modes(Kr, Mr, 4);
check(sprintf('Rigid surface (%d dofs) == meshed rigid body', size(Kr, 1)), f, fRig, 1e-8);
a = fembd.Assembly();  a.add(L);
a.add(fembd.Body.pointMass(mR, [3*Lx/4 h/2 h/2], JR, Name = 'tipmass'));
r = a.couple('left', @(X) X(:,1) > Lx/2 - 1e-9);
a.joint(r, 'tipmass', 'fixed');  a.fix('left', @(X) X(:,1) < 1e-9);
[Kr, Mr] = a.system();  [~, f] = fembd.modes(Kr, Mr, 4);
check('Point mass + fixed joint == meshed rigid body', f, fRig, 1e-6);

% joints: two rigid quarters joined by a 'fixed' joint == one rigid half
fembd.boxSTL(fullfile(tmp, 'q1.stl'), [Lx/2 0 0],   [3*Lx/4 h h]);
fembd.boxSTL(fullfile(tmp, 'q2.stl'), [3*Lx/4 0 0], [Lx h h]);
Q1 = fembd.Body(fullfile(tmp, 'q1.stl'), steel, Size = h, Order = 1, Rigid = true);
Q2 = fembd.Body(fullfile(tmp, 'q2.stl'), steel, Size = h, Order = 1, Rigid = true);
a = fembd.Assembly();  a.add(L);  a.add(Q1);  a.add(Q2);
a.tie('left', 'q1');
a.joint('q1', 'q2', 'fixed', Point = [3*Lx/4 h/2 h/2]);
a.fix('left', @(X) X(:,1) < 1e-9);
[Kr, Mr] = a.system();  [~, f] = fembd.modes(Kr, Mr, 4);
check('Fixed joint of two rigid parts == one rigid part', f, fRig, 1e-8);

% revolute and prismatic joints to ground with springs: 1-dof oscillators
I = m*(Lx^2 + h^2)/12 + m*(Lx/2)^2;                    % about a z-axis at x = 0
a = fembd.Assembly();  a.add(Rb);
a.joint('ground', 'beam', 'revolute', Point = [0 h/2 h/2], Axis = [0 0 1], Stiffness = 1e3);
[Kr, Mr] = a.system();
check('Revolute joint + torsion spring: sqrt(k/I)', [size(Kr, 1), sqrt(Kr/Mr)], [1, sqrt(1e3/I)], 1e-6);
a = fembd.Assembly();  a.add(Rb);
a.joint('ground', 'beam', 'prismatic', Axis = [1 1 0], Stiffness = 1e6);
[Kr, Mr] = a.system();
check('Prismatic joint + spring: sqrt(k/m)', [size(Kr, 1), sqrt(Kr/Mr)], [1, sqrt(1e6/m)], 1e-6);

% coupling (RBE2) of the clamped face, grounded by a fixed joint == fix
a = fembd.Assembly();  a.add(L);
r = a.couple('left', @(X) X(:,1) < 1e-9);
a.joint('ground', r, 'fixed');
[Kr, Mr] = a.system();  [~, f1] = fembd.modes(Kr, Mr, 4);
a = fembd.Assembly();  a.add(L);  a.fix('left', @(X) X(:,1) < 1e-9);
[Kr, Mr] = a.system();  [~, f2] = fembd.modes(Kr, Mr, 4);
check('Coupling + fixed joint to ground == clamp', f1, f2, 1e-10);

% revolute joint on a coupled elastic end: pinned-free beam, one zero mode
a = fembd.Assembly();  a.add(L);
r = a.couple('left', @(X) X(:,1) < 1e-9);
j = a.joint('ground', r, 'revolute', Axis = [0 1 0]);
[Kr, Mr] = a.system();  [~, f] = fembd.modes(Kr, Mr, 3);
beam = @(lam) lam^2/(2*pi)*sqrt(210e9*h^4/12/(7800*h^2*(Lx/2)^4));
% about y the beam is pinned-free, about z still clamped-free; at L/h = 10
% shear and rotary inertia lower the pinned-free mode ~2.4% below Euler-Bernoulli
check('Hinged beam: 0 Hz, clamped-free (z), pinned-free (y)', ...
      [f(1)/f(2), f(2:3)'], [0, beam(1.875104), beam(3.926602)], 3e-2);
check('Joint coordinate is a dof of q', numel(a.qdof(a.dof(j))), 1, 0);

% posing joints: rigid two-link arm at a 60 degree elbow vs the analytic
% double-pendulum mass matrix (coordinates: shoulder angle, elbow angle)
fembd.boxSTL(fullfile(tmp, 'l1.stl'), [0 -0.01 -0.005],   [0.2 0.01 0.005]);
fembd.boxSTL(fullfile(tmp, 'l2.stl'), [0.2 -0.01 -0.005], [0.4 0.01 0.005]);
a = fembd.Assembly();
a.add(fembd.Body.rigidSurface(fullfile(tmp, 'l1.stl'), 2700));
a.add(fembd.Body.rigidSurface(fullfile(tmp, 'l2.stl'), 2700));
kj = 200;
a.joint('ground', 'l1', 'revolute', Point = [0 0 0],   Stiffness = kj, Name = 'shoulder');
a.joint('l1', 'l2',     'revolute', Point = [0.2 0 0], Stiffness = kj, Name = 'elbow');
X0 = a.X;  phi = pi/3;
a.pose('elbow', phi);
[Kr, Mr] = a.system();  [~, f] = fembd.modes(Kr, Mr, 2);
mL = 2700*0.2*0.02*0.01;  Ic = mL*(0.2^2 + 0.02^2)/12;  I2 = Ic + mL*0.1^2;
Mq = [Ic + mL*0.1^2 + I2 + mL*0.2^2 + 2*mL*0.2*0.1*cos(phi), I2 + mL*0.2*0.1*cos(phi);
      I2 + mL*0.2*0.1*cos(phi), I2];
check('Posed rigid arm (60 deg elbow) vs analytic', f, sqrt(sort(eig(kj*eye(2), Mq)))/(2*pi), 1e-8);
a.pose('shoulder', 0.7);
[Kr, Mr] = a.system();  [~, f2] = fembd.modes(Kr, Mr, 2);
check('Posing the shoulder: frequencies invariant', f2, f, 1e-10);
a.pose('shoulder', 0);  a.pose('elbow', 0);
check('pose(joint, 0) restores the geometry', max(abs(a.X(:) - X0(:))), 0, 1e-14);

% prismatic pose carries the downstream joint along
a = fembd.Assembly();
a.add(fembd.Body.rigidSurface(fullfile(tmp, 'l1.stl'), 2700, Name = 'slider'));
a.add(fembd.Body.pointMass(1, [0.3 0 0], 1e-3*eye(3)));
a.joint('ground', 'slider', 'prismatic', Axis = [1 0 0], Name = 'rail');
a.joint('slider', 'mass', 'revolute', Point = [0.3 0 0], Name = 'pin');
a.pose('rail', 0.05);
check('Prismatic pose moves bodies and child joints', ...
      [a.Bodies{2}.Ref(1), a.Joints(2).p(1), a.Bodies{1}.Ref(1)], [0.35 0.35 0.15], 1e-7);   % STL is single precision

% a closed loop (door on two hinges) cannot be posed
a = fembd.Assembly();  a.add(L);
r1 = a.couple('left', @(X) X(:,1) < 1e-9 & X(:,3) < h/2);
r2 = a.couple('left', @(X) X(:,1) < 1e-9 & X(:,3) > h/2);
a.joint('ground', r1, 'revolute', Axis = [0 0 1]);
a.joint('ground', r2, 'revolute', Axis = [0 0 1]);
try, a.pose('joint1', 0.1); ok = false; catch e, ok = strcmp(e.identifier, 'fembd:pose'); end
check('Posing a closed loop gives an error', ok, true, 0);

% closed loops: planar four-bar of point-mass links (revolute axes along z,
% so the loop is also overconstrained in 3D). Reference: velocity analysis.
O2 = [0 0];  A = [0.05 0.1];  B = [0.3 0.15];  O4 = [0.25 0];
mk = [0.2 0.5 0.3];  Jz = [1e-4 4e-4 2e-4];  kj = 50;
cs = {(O2 + A)/2, (A + B)/2, (B + O4)/2};
perp = @(v) [-v(2), v(1)];
vA = perp(A - O2);                                    % crank at 1 rad/s
w34 = [perp(B - A)', -perp(B - O4)'] \ (-vA');        % [w3; w4]
vc3 = vA + w34(1)*perp(cs{2} - A);
Ieff = Jz(1) + mk(1)*norm(cs{1} - O2)^2 + mk(2)*norm(vc3)^2 + Jz(2)*w34(1)^2 ...
     + (Jz(3) + mk(3)*norm(cs{3} - O4)^2)*w34(2)^2;
fourbar = @(order, kSpring) buildFourbar(order, kSpring, mk, Jz, cs, O2, A, B, O4);
[a, Kr, Mr] = fourbar([1 2 3 4], [kj 0 0 0]);
check('Four-bar: 1 dof, spring on crank vs velocity analysis', ...
      [size(Kr, 1), full(sqrt(Kr/Mr))/(2*pi)], [1, sqrt(kj/Ieff)/(2*pi)], 1e-8);
check('Four-bar: crank coordinate stays independent', numel(a.qdof(a.dof('j1'))), 1, 0);
[~, Kr, Mr] = fourbar([1 2 3 4], [0 0 0 kj]);         % spring in the cut joint
check('Four-bar: spring on the cut joint (k w4^2)', full(sqrt(Kr/Mr))/(2*pi), ...
      sqrt(kj*w34(2)^2/Ieff)/(2*pi), 1e-8);
[~, Kr, Mr] = fourbar([4 3 2 1], [kj 0 0 0]);         % other joint is cut
check('Four-bar: independent of joint order', full(sqrt(Kr/Mr))/(2*pi), sqrt(kj/Ieff)/(2*pi), 1e-8);

% a bar hinged to ground at both ends is locked: no dofs remain
a = fembd.Assembly();  a.add(fembd.Body.pointMass(1, [0.5 0 0], eye(3)*1e-2, Name = 'bar'));
a.joint('ground', 'bar', 'revolute', Point = [0 0 0]);
a.joint('ground', 'bar', 'revolute', Point = [1 0 0]);
[Kr, ~] = a.system();
check('Bar hinged at both ends: no dofs', size(Kr, 1), 0, 0);

% transmissions: motor rotor -> rz of a spherical joint, ratio N
mL = 3;  cL = [0.17 0 0];  JL = 0.02*eye(3);  Jm = 1e-4;  N = 2;  kb = 169430;  kz = 500;
% (centre of mass in the xy plane: rz is decoupled from the tilt directions)
Js = JL(3,3) + mL*(cL(1)^2 + cL(2)^2);
[~, Kr, Mr] = buildBelt(kb, N, mL, cL, JL, Jm, 0);              % elastic belt
[~, f] = fembd.modes(Kr, Mr, size(Kr, 1));                       % rigid, 2 tilt, belt mode
check('Elastic transmission vs two-inertia drive mode', f(end), ...
      sqrt(kb*(1/Js + 1/(N^2*Jm)))/(2*pi), 1e-8);
[a, Kr, Mr] = buildBelt([], N, mL, cL, JL, Jm, kz);            % rigid gear
[~, f] = fembd.modes(Kr, Mr, 1);
check('Rigid transmission: reflected inertia Js + N^2 Jm', f, sqrt(kz/(Js + N^2*Jm))/(2*pi), 1e-8);
check('Rigid transmission: driver coordinate stays in q', numel(a.qdof(a.dof('motor'))), 1, 0);

% transmissions on reference points (absolute rotations of rigid bodies)
f1 = refBelt(1, 6, kb);  f2 = refBelt(0, 1, kb);      % link rz vs its joint coordinate
check('Transmission on a rigid body rz == on its joint', f1, f2, 1e-10);
f1 = refBelt(1, 6, []);  f2 = refBelt(0, 1, []);
check('Rigid transmission on a rigid body rz == on its joint', f1, f2, 1e-10);
f1 = elbowBelt(true, kb);  f2 = elbowBelt(false, kb);  % base-mounted elbow drive
check('Base-mounted elbow belt: absolute rz == shoulder + elbow', f1, f2, 1e-10);

% setMass: mass properties of rigid bodies from CAD data
B = fembd.Body.rigidSurface(fullfile(tmp, 'beam.stl'), 7800);
c0 = B.MassProps.c;  J0 = B.MassProps.J;  m0 = B.MassProps.m;
B.setMass(2*m0);
check('setMass(m): mass set, inertia scaled, c kept', [B.MassProps.m, B.MassProps.J(:)', B.MassProps.c], ...
      [2*m0, 2*J0(:)', c0], 1e-12);
Jo = J0 + m*(dot(c0, c0)*eye(3) - c0'*c0);            % the same inertia about the origin
B.setMass(m, CenterOfMass = c0, Inertia = Jo, InertiaPoint = [0 0 0]);
check('setMass: inertia about another point (parallel axis)', B.MassProps.J, J0, 1e-10);
Rz = [cosd(30) -sind(30) 0; sind(30) cosd(30) 0; 0 0 1];
B.setMass(m, Inertia = [1 2 2.5], Axes = Rz);
check('setMass: principal moments and axes', B.MassProps.J, Rz*diag([1 2 2.5])*Rz', 1e-12);
B.setMass(m*1, CenterOfMass = 1e3*c0, Inertia = 1e6*J0, LengthScale = 1e-3);
check('setMass: mm and kg mm^2', [B.MassProps.c, B.MassProps.J(:)'], [c0, J0(:)'], 1e-12);
try, B.setMass(m, Inertia = [1 1 3]); ok = false; catch e, ok = strcmp(e.identifier, 'fembd:setMass'); end
check('setMass: invalid inertia (triangle inequality) is rejected', ok, true, 0);
a = fembd.Assembly();  a.add(Rb);  Rb.setMass(3*m);       % meshed rigid body: no double counting
[~, Mr] = a.system();
check('setMass on a meshed rigid body replaces the mesh mass', full(Mr(1,1)), 3*m, 1e-10);
Rb.MassProps = [];                                       % restore (Rb is used elsewhere)

% damping on joints and transmissions
mD = 2;  rD = [0.2 0 0];  JD = 0.01;  kD = 300;  cD = 0.4;
a = fembd.Assembly();
a.add(fembd.Body.pointMass(mD, rD, JD*eye(3), Name = 'link'));
a.joint('ground', 'link', 'revolute', Point = [0 0 0], Stiffness = kD, Damping = cD, Name = 'axis');
a.input('tau', 'axis');  a.output('th', 'axis');
[~, ~, ~, Dr] = a.system();
check('Joint damping: D = c on the joint coordinate', full(Dr), cD, 1e-12);
sys = a.statespace(Modes = 1, Damping = 0);
Jtot = JD + mD*norm(rD)^2;
check('Joint damping: poles of J s^2 + c s + k', sort(eig(sys.A)), sort(roots([Jtot cD kD])), 1e-10);

a = fembd.Assembly();                                     % damped belt between rotor and link
a.add(fembd.Body.pointMass(2, [0.2 0 0], 0.01*eye(3), Name = 'link'));
a.add(fembd.Body.pointMass(0.1, [0.5 0 0], 1e-4*eye(3), Name = 'rotor'));
a.joint('ground', 'link',  'revolute', Point = [0 0 0], Stiffness = 300, Name = 'axis');
a.joint('ground', 'rotor', 'revolute', Point = [0.5 0 0], Damping = 1e-3, Name = 'motor');
a.transmission({'motor', 'link'}, [1/3 -1], Stiffness = 2e4, Damping = 2, Dirs = {1, 6});
a.input('tau', 'motor');  a.output('th', 'link', [], 6);
[Kr, Mr, ~, Dr] = a.system();
sys = a.statespace(Modes = size(Kr, 1), Damping = 0);
pe = polyeig(full(Kr), full(Dr), full(Mr));
check('Transmission + joint damping: statespace poles == polyeig(K, D, M)', ...
      sort(eig(sys.A), 'ComparisonMethod', 'real'), sort(pe, 'ComparisonMethod', 'real'), 1e-8);

% plot handles and mode animation (headless): beam + rigid block + point mass
a = fembd.Assembly();  a.add(L);  a.add(fembd.Body.rigidSurface(fullfile(tmp, 'right.stl'), 7800));
a.add(fembd.Body.pointMass(0.5, [Lx/4 h h], 1e-4*eye(3), Name = 'pm'));
a.tie('left', 'right');  a.fix('left', @(X) X(:,1) < 1e-9);
r = a.couple('left', @(X) abs(X(:,1) - Lx/4) < h/4);  a.joint(r, 'pm', 'fixed');
[Kr, Mr, T] = a.system();  [V, f] = fembd.modes(Kr, Mr, 1);
fig = figure('Visible', 'off');
hp = a.plot(T*V(:,1), 0.01);
check('plot returns one handle per body (2 patches, 1 marker)', ...
      [numel(hp), sum(isa(hp(1), 'matlab.graphics.primitive.Patch') + isa(hp(2), 'matlab.graphics.primitive.Patch')), isa(hp(3), 'matlab.graphics.chart.primitive.Line')], [3 2 1], 0);
tic;  a.animate(T*V(:,1), Frequency = f(1), Rate = 2, Duration = 0.4);  ta = toc;
hp = findobj(gca, 'Type', 'patch');  hm = findobj(gca, 'Type', 'line');
moved = any(abs(hp(end).Vertices(:) - L.X(:)) > 0) && ~isequal([hm.XData hm.YData hm.ZData], [Lx/4 h h]);
check(sprintf('animate: runs for its Duration (%.2f s) and moves bodies and markers', ta), ...
      double(ta >= 0.4 && ta < 5 && moved && contains(get(get(gca, 'Title'), 'String'), 'slow motion')), 1, 0);
close(fig);

% memory: an elastic transmission on a model with ~600k dofs (a detailed rigid
% STL) must not build full-length outer products (g'*g ran out of memory)
fembd.cylinderSTL(fullfile(tmp, 'big.stl'), [0 0], 0.03, 0, 0.3, 100000);
a = fembd.Assembly();
a.add(fembd.Body.rigidSurface(fullfile(tmp, 'big.stl'), 7850, Name = 'link'));
a.add(fembd.Body.pointMass(0.1, [0.1 0 0], 1e-4*eye(3), Name = 'rotor'));
a.joint('ground', 'link',  'revolute', Point = [0 0 0],   Name = 'axis');
a.joint('ground', 'rotor', 'revolute', Point = [0.1 0 0], Name = 'motor');
a.transmission({'motor', 'link'}, [1 -1], Stiffness = 169430, Dirs = {1, 6});
tic;  [Kr, ~] = a.system();  t = toc;
check(sprintf('Transmission at %d dofs: system in %.1f s', a.nDof, t), size(Kr, 1), 2, 0);

% inputs, outputs and state space: cantilever, tip force -> tip deflection
a = fembd.Assembly();  a.add(L);  a.fix('left', @(X) X(:,1) < 1e-9);
a.input('Fz', 'left', @(X) X(:,1) > Lx/2 - 1e-9, 3);
a.output('uz', 'left', @(X) X(:,1) > Lx/2 - 1e-9, 3);
[Kr, Mr, T] = a.system();  [Bu, Cu] = a.ioMatrices();
dc = Cu*T*(Kr\(T'*Bu));                                         % exact static gain
sys = a.statespace(Modes = 2);
dcSys = -sys.C*(sys.A\sys.B) + sys.D;
check('State space: DC gain with static correction (2 modes)', dcSys, dc, 1e-8);
sys0 = a.statespace(Modes = 2, StaticCorrection = false);
err0 = abs(-sys0.C*(sys0.A\sys0.B)/dc - 1);
check(sprintf('State space: correction matters (error %.1e without)', err0), err0 > 1e-4, true, 0);
[~, f] = fembd.modes(Kr, Mr, 2);
check('State space: modal frequencies', sys.f(1:2), f, 1e-8);

% mechanism (free drive) and acceleration feedthrough: belt drive, all modes
[a, Kr, Mr] = buildBelt(kb, N, mL, cL, JL, Jm, 0);
a.input('tau', 'motor');
a.output('acc', 'motor', [], [], Type = 'acceleration');
a.output('th', 'shoulder', [], 3);
[~, ~, T] = a.system();  [Bu, Cu] = a.ioMatrices();
sys = a.statespace(Modes = size(Kr, 1), Damping = 0);
check('State space: free drive gives a double integrator', sum(abs(eig(sys.A)) < 1e-6), 2, 0);
check('State space: acceleration feedthrough C M^-1 B', sys.D(1), full(Cu(1,:)*T*(Mr\(T'*Bu))), 1e-8);
check('State space: signal names', [sys.InputName, sys.OutputName], {'tau', 'acc', 'th'}, 0);

% MSH 4.1 reader vs 2.2
geo = fullfile(tmp, 'b.geo');
writelines(sprintf(['Merge "%s";\nClassifySurfaces{40*Pi/180,1,1,Pi};\nCreateGeometry;\n' ...
    'Surface Loop(1) = Surface{:};\nVolume(1) = {1};\n'], strrep(fullfile(tmp, 'beam.stl'), '\', '/')), geo);
for v = ["22" "41"]
    fembd.gmsh(sprintf('"%s" -3 -order 2 -clmax %g -format msh%s -o "%s/b%s.msh" -v 0', ...
        geo, h, v, tmp, v));
end
[X2, T2] = fembd.readMsh(fullfile(tmp, 'b22.msh'));
[X4, T4] = fembd.readMsh(fullfile(tmp, 'b41.msh'));
check('MSH 4.1 == MSH 2.2', norm(X2(T2(:)) - X4(T4(:))) + abs(numel(T2) - numel(T4)), 0, 1e-12);

% nonlinear model: double compound pendulum vs Lagrange equations (ode45)
m1 = 1.5;  m2 = 0.8;  L1 = 0.4;  lc1 = 0.2;  lc2 = 0.25;  I1 = 0.02;  I2 = 0.01;  g = 9.81;
a = fembd.Assembly();
a.add(fembd.Body.pointMass(m1, [lc1 0 0], I1*eye(3), Name = 'link1'));
a.add(fembd.Body.pointMass(m2, [L1+lc2 0 0], I2*eye(3), Name = 'link2'));
a.joint('ground', 'link1', 'revolute', Point = [0 0 0], Axis = [0 1 0], Name = 'j1');
a.joint('link1', 'link2', 'revolute', Point = [L1 0 0], Axis = [0 1 0], Name = 'j2');
a.output('ry1', 'link1', [], 5);  a.output('ry2', 'link2', [], 5);  a.output('q2', 'j2');
[t, y] = a.nonlinear(Gravity = [0 0 -g]).simulate([0 1], Step = 2.5e-4);
c12 = @(z) m2*L1*lc2*cos(z(1) - z(2));  s12 = @(z) m2*L1*lc2*sin(z(1) - z(2));
f = @(t, z) [z(3:4); [I1 + m1*lc1^2 + m2*L1^2, c12(z); c12(z), I2 + m2*lc2^2] \ ...
    -[s12(z)*z(4)^2 + (m1*lc1 + m2*L1)*g*sin(z(1)); -s12(z)*z(3)^2 + m2*lc2*g*sin(z(2))]];
[~, z] = ode45(f, t, [pi/2; pi/2; 0; 0], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
check('Nonlinear: double pendulum vs Lagrange (ode45)', ...
      pi/2 - unwrap([y.ry1 y.ry2]), z(:,1:2), 2e-5);
check('Nonlinear: joint angle == difference of link angles', wrapToPi(y.q2 - (y.ry2 - y.ry1)), 0, 1e-12);

% nonlinear model: damped elastic belt at large angles == linear (parallel axes)
[a0, ~] = refBeltNL(false);  [a, ~] = refBeltNL(true);
[K, M, ~, D] = a0.system();  n = size(K, 1);
i = [a0.qdof(a0.dof('motor')), a0.qdof(a0.dof('axis'))];
A = [zeros(n), eye(n); -M\K, -M\D];  x0 = zeros(2*n, 1);  x0(i) = [1.5; 0.5];
[t, y] = a.nonlinear().simulate([0 0.3], Step = 1e-4);
X = cell2mat(arrayfun(@(tk) expm(A*tk)*x0, t', 'UniformOutput', false))';
check('Nonlinear: belt, 1.5 rad, == linear solution', [y.qm y.qa], X(:, i), 2e-4);

% nonlinear model: rigid transmission, motor torque, redundant support
N = 3;  Jm = 1e-4;  Js = 0.01 + 2*0.2^2;  tau = 0.02;
a = fembd.Assembly();
a.add(fembd.Body.pointMass(2, [0.2 0 0], 0.01*eye(3), Name = 'link'));
a.add(fembd.Body.pointMass(0.1, [0.5 0 0], Jm*eye(3), Name = 'rotor'));
a.joint('ground', 'link',  'revolute', Point = [0 0 0], Name = 'axis');
a.joint('ground', 'rotor', 'revolute', Point = [0.5 0 0], Name = 'motor');
a.transmission({'motor', 'axis'}, [1/N -1]);
a.fix('rotor', [], 1:3);                                 % redundant with the motor joint
a.input('tau', 'motor');  a.output('qm', 'motor');  a.output('qa', 'axis');
nl = a.nonlinear();
[t, y] = nl.simulate([0 0.5], Step = 1e-3, Input = @(t, y) tau);
check('Nonlinear: rigid transmission + torque, (Jm + Js/N^2) qm'''' = tau', ...
      [y.qm, N*y.qa], repmat(0.5*tau/(Jm + Js/N^2)*t.^2, 1, 2), 1e-7);
check('Nonlinear: redundant constraints dropped (14 -> 11 rows)', numel(nl.keep), 11, 0);

% nonlinear model: heavy top on a spherical joint (energy, ode45)
a = fembd.Assembly();
a.add(fembd.Body.pointMass(1, [0 0 0.1], diag([2 2 4])*1e-3, Name = 'top'));
a.joint('ground', 'top', 'spherical', Point = [0 0 0], Name = 'pivot');
a.pose('pivot', [0.3 0 0]);
nl = a.nonlinear(Gravity = [0 0 -g]);
[t, ~, x] = nl.simulate([0 0.2], Step = 1e-4, Rates = struct('pivot', [0 0 60]), Rho = 0.9);
Jb = nl.Bodies(1).J;  w = x.nu(4:6,:);  r = squeeze(x.r(:,1,:));
E = 0.5*sum(x.nu(1:3,:).^2) + 0.5*sum(w.*(Jb*w)) + g*r(3,:);
cb = r(:,1);  Jp = Jb + (cb'*cb*eye(3) - cb*cb');      % body frame = world at t = 0
f = @(t, z) [reshape(reshape(z(1:9), 3, 3)*[0 -z(12) z(11); z(12) 0 -z(10); -z(11) z(10) 0], 9, 1); ...
             Jp\(cross(cb, reshape(z(1:9), 3, 3)'*[0; 0; -g]) - cross(z(10:12), Jp*z(10:12)))];
[~, z] = ode45(f, [0 0.2], [reshape(eye(3), 9, 1); w(:,1)], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
check(sprintf('Nonlinear: heavy top (%.0f rad/s), energy', norm(w(:,1))), E/E(1), 1, 1e-6);
check('Nonlinear: heavy top vs Euler equations (ode45)', r(:,end), reshape(z(end,1:9), 3, 3)*cb, 1e-6);

% nonlinear model: free fall along a prismatic joint
a = fembd.Assembly();
a.add(fembd.Body.pointMass(3, [0.1 0.2 0.3], 0.01*eye(3), Name = 'slider'));
a.joint('ground', 'slider', 'prismatic', Axis = [0 0 1], Name = 'z');
a.output('z', 'z');
[t, y] = a.nonlinear(Gravity = [1 2 -g]).simulate([0 0.3], Step = 1e-2);
check('Nonlinear: free fall on a prismatic joint', y.z, -0.5*g*t.^2, 1e-12);

% nonlinear elastic: cantilever with a tip mass on its second coupling == linear
Lb = 0.5;  hb = 0.02;  beam = fullfile(tmp, 'nlbeam.stl');
fembd.boxSTL(beam, [0 0 0], [Lb hb hb]);
a = elasticBeam(beam, steel, Lb, [0 2e-5]);
a.add(fembd.Body.pointMass(0.3, [Lb + 0.03, hb/2, hb/2], 1e-4*eye(3), Name = 'mass'));
a.joint('ground', 'root', 'fixed');  a.joint('tip', 'mass', 'fixed', Point = [Lb hb/2 hb/2]);
a.input('Fz', 'mass', [], 3);  a.output('uz', 'mass', [], 3);
a.output('mid', 'beam', @(X) abs(X(:,1) - Lb/2) < 0.006, 3);
sys = a.statespace(Modes = 30, Damping = 0);
[t, y] = a.nonlinear(Modes = 10).simulate([0 0.02], Step = 4e-5, Input = @(t, y) 1, Rho = 0.95);
Ad = expm(sys.A*(t(2) - t(1)));  Bd = sys.A\((Ad - eye(size(Ad)))*sys.B);
xs = zeros(size(sys.A, 1), 1);  Yl = zeros(numel(t), 2);
for k = 2:numel(t), xs = Ad*xs + Bd;  Yl(k,:) = (sys.C*xs)'; end
check('Nonlinear elastic: tip mass and node set == linear (step force)', ...
      [y.uz, y.mid]./max(abs(Yl)), Yl./max(abs(Yl)), 1e-3);

% nonlinear elastic: stiff flexible pendulum == rigid pendulum (3 rad swing)
stiff = fembd.material(210e12, 0.3, 7800);
[~, yr] = elasticPendulum(beam, stiff, Lb, true).nonlinear(Gravity = [0 0 -g]).simulate([0 0.6]);
[~, ye] = elasticPendulum(beam, stiff, Lb, false).nonlinear(Gravity = [0 0 -g], Modes = 6).simulate([0 0.6]);
check(sprintf('Nonlinear elastic: stiff pendulum (%.1f rad) == rigid', max(abs(yr.th))), ye.th, yr.th, 1e-4);

% nonlinear elastic: soft flexible pendulum conserves energy; kinetic energy
% from the inertia invariants == nodal 1/2 xdot' M xdot
nl = elasticPendulum(beam, fembd.material(2e9, 0.3, 7800), Lb, false).nonlinear(Gravity = [0 0 -g], Modes = 6);
[~, ~, x] = nl.simulate([0 0.6], Step = 5e-4, Rho = 1);
E = nl.energy(x);  B = nl.Bodies(1);  k = size(x.nu, 2);
U = B.Xl' + reshape(B.Phi*x.e(:,k), 3, []);  w = x.nu(4:6,k);
xd = x.nu(1:3,k) + x.R(:,:,1,k)*(cross(repmat(w, 1, size(U, 2)), U) + reshape(B.Phi*x.nu(7:end,k), 3, []));
[~, Mb] = nl.Assembly.Bodies{1}.matrices();
check(sprintf('Nonlinear elastic: soft pendulum energy (strain %.0e J)', max(E.strain)), ...
      (E.total - E.total(1))/max(E.kinetic), 0, 1e-5);
check('Nonlinear elastic: kinetic energy, invariants == nodal', E.kinetic(k), 0.5*full(xd(:)'*Mb*xd(:)), 1e-10);
end

function a = elasticBeam(file, mat, L, damping)
% Elastic Tet4 beam along x with couplings at both ends.
a = fembd.Assembly();
a.add(fembd.Body(file, mat, Size = 0.01, Order = 1, Name = 'beam', Damping = damping));
a.couple('beam', @(X) X(:,1) < 1e-9, Name = 'root');
a.couple('beam', @(X) X(:,1) > L - 1e-9, Name = 'tip');
end

function a = elasticPendulum(file, mat, L, rigid)
% Beam pendulum hinged about y at its root, rigid or elastic.
if rigid
    a = fembd.Assembly();
    a.add(fembd.Body(file, mat, Size = 0.01, Order = 1, Name = 'beam', Rigid = true));
    a.joint('ground', 'beam', 'revolute', Point = [0 0.01 0.01], Axis = [0 1 0], Name = 'hinge');
else
    a = elasticBeam(file, mat, L, [0 0]);
    a.joint('ground', 'root', 'revolute', Point = [0 0.01 0.01], Axis = [0 1 0], Name = 'hinge');
end
a.output('th', 'hinge');
end

function [a, name] = refBeltNL(posed)
% Link on a torsion spring driven by a rotor through a damped elastic belt.
a = fembd.Assembly();  name = 'belt';
a.add(fembd.Body.pointMass(2, [0.2 0 0], 0.01*eye(3), Name = 'link'));
a.add(fembd.Body.pointMass(0.1, [0.5 0 0], 1e-4*eye(3), Name = 'rotor'));
a.joint('ground', 'link',  'revolute', Point = [0 0 0], Stiffness = 300, Damping = 0.05, Name = 'axis');
a.joint('ground', 'rotor', 'revolute', Point = [0.5 0 0], Damping = 1e-3, Name = 'motor');
a.transmission({'motor', 'axis'}, [1/3 -1], Stiffness = 50, Damping = 0.02, Name = name);
a.output('qm', 'motor');  a.output('qa', 'axis');
if posed, a.pose('motor', 1.5);  a.pose('axis', 0.5); end
end

function [a, K, M] = buildFourbar(order, k, mk, Jz, cs, O2, A, B, O4)
% Four-bar of point masses; joints added in the given order.
a = fembd.Assembly();
for i = 1:3
    a.add(fembd.Body.pointMass(mk(i), [cs{i} 0], diag([1 1 1]*Jz(i)), Name = sprintf('link%d', i)));
end
J = {{'ground', 'link1', O2}, {'link1', 'link2', A}, {'link2', 'link3', B}, {'ground', 'link3', O4}};
for i = order
    a.joint(J{i}{1}, J{i}{2}, 'revolute', Point = [J{i}{3} 0], Stiffness = k(i), ...
            Name = sprintf('j%d', i));
end
[K, M] = a.system();
end

function [a, K, M] = buildBelt(k, N, mL, cL, JL, Jm, kz)
% Link on a spherical joint (stiff tilt, rz spring kz) driven by a motor rotor
% on a revolute joint through a transmission theta_m/N - theta_s.
a = fembd.Assembly();
a.add(fembd.Body.pointMass(mL, cL, JL, Name = 'link'));
a.add(fembd.Body.pointMass(0.1, [0.1 0 0], Jm*eye(3), Name = 'rotor'));
a.joint('ground', 'link', 'spherical', Point = [0 0 0], Stiffness = [1e5 1e5 kz], Name = 'shoulder');
a.joint('ground', 'rotor', 'revolute', Point = [0.1 0 0], Axis = [0 0 1], Name = 'motor');
a.transmission({'motor', 'shoulder'}, [1/N -1], Stiffness = k, Dirs = {1, 3});
[K, M] = a.system();
end

function f = refBelt(onBody, dir, k)
% Motor driving a link on a revolute joint, via the link's rz or the joint.
a = fembd.Assembly();
a.add(fembd.Body.pointMass(2, [0.2 0 0], 0.01*eye(3), Name = 'link'));
a.add(fembd.Body.pointMass(0.1, [0.5 0 0], 1e-4*eye(3), Name = 'rotor'));
a.joint('ground', 'link',  'revolute', Point = [0 0 0], Stiffness = 300, Name = 'axis');
a.joint('ground', 'rotor', 'revolute', Point = [0.5 0 0], Name = 'motor');
ids = {'motor', 'axis'};  if onBody, ids{2} = 'link'; end
a.transmission(ids, [1/3 -1], Stiffness = k, Dirs = {1, dir});
[K, M] = a.system();  [~, f] = fembd.modes(K, M, size(K, 1));
end

function f = elbowBelt(absolute, k)
% Two-link arm (spherical joints, free rz) with the elbow driven from the
% base: via the lower link's absolute rz, or via shoulder + elbow coordinates.
a = fembd.Assembly();
a.add(fembd.Body.pointMass(3, [0.17 0 0], 0.02*eye(3), Name = 'upper'));
a.add(fembd.Body.pointMass(2, [0.5 0 0],  0.01*eye(3), Name = 'lower'));
a.add(fembd.Body.pointMass(0.1, [0 0.1 0], 2e-4*eye(3), Name = 'rotor'));
tilt = [1e5 1e5 0];
a.joint('ground', 'upper', 'spherical', Point = [0 0 0],    Stiffness = tilt + [0 0 400], Name = 'shoulder');
a.joint('upper',  'lower', 'spherical', Point = [0.34 0 0], Stiffness = tilt, Name = 'elbow');
a.joint('ground', 'rotor', 'revolute',  Point = [0 0.1 0],  Name = 'motor2');
if absolute
    a.transmission({'motor2', 'lower'}, [1/4 -1], Stiffness = k, Dirs = {1, 6});
else
    a.transmission({'motor2', 'shoulder', 'elbow'}, [1/4 -1 -1], Stiffness = k, Dirs = {1, 3, 3});
end
[K, M] = a.system();  [~, f] = fembd.modes(K, M, size(K, 1));
end

function check(name, val, ref, tol)
if iscell(val), val = double(isequal(val, ref)); ref = 1; end
err = max(abs(val(:) - ref(:))./max(abs(ref(:)), 1));
if isequal(ref, 0), err = max(abs(val(:))); end
status = "PASS";  if err > tol, status = "FAIL"; end
fprintf('%-4s  %-45s err = %.2e\n', status, name, err);
end
