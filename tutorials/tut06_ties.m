%% Tutorial 6 - Ties, elastic interfaces and supports
% Each body is meshed separately, so the meshes do not match at an
% interface. This tutorial covers:
%
% * bonded ties between non-matching meshes, and elastic groups
% * choosing master and slave, and the |Tol| and |Select| options
% * elastic interfaces (distributed springs) with |Stiffness|
% * ties between elastic and rigid bodies
% * clamping elastic bodies on nodes

addpath(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(tempdir, 'fmbdTutorials');  if ~isfolder(out), mkdir(out); end
steel = fmbd.material(210e9, 0.3, 7850);
clamp = @(X) X(:,1) < 1e-6;
fmbd.boxSTL(fullfile(out, 'beam.stl'),  [0 0 0],   [1000 50 50]);
fmbd.boxSTL(fullfile(out, 'left.stl'),  [0 0 0],   [500 50 50]);
fmbd.boxSTL(fullfile(out, 'right.stl'), [500 0 0], [1000 50 50]);
body = @(name, h) fmbd.Body.elastic(fullfile(out, [name '.stl']), steel, Size = h, ...
                                    Scale = 1e-3, Name = name, Modes = 10);

%% 1. Clamping on nodes
% |s.fix(body, sel, dirs)| clamps the selected nodes of an elastic body (all
% three directions by default). The body is then grounded: its frame is
% fixed, and its modes are those of the supported structure. The supports
% must hold the body (e.g. a clamped face); for a sliding or rotating
% support use a coupling and a joint. A one-piece cantilever as reference:
s = fmbd.System();
s.add(body('beam', 25));
s.fix('beam', clamp);
fRef = s.modes([], 6);

%% 2. Bonded tie between non-matching meshes
% The same beam as two halves with different element sizes.
% |s.tie(master, slave)| takes the slave's surface nodes that lie within Tol
% of the master surface, and constrains each one to the point of the master
% face it touches: u_slave = sum_j N_j(xi) u_master,j. Use the finer mesh as
% slave, so that every slave node has a master face to follow. Tied elastic
% bodies become one elastic group: one frame, and modes of the combined
% structure (the Modes of the bodies add up).
s = fmbd.System();
s.add(body('left', 25));
s.add(body('right', 15));
ns = s.tie('left', 'right');
s.fix('left', clamp);
[f, V] = s.modes([], 6);
fprintf('%d slave nodes tied\n', ns);
disp(table(fRef, f, 'VariableNames', {'one_piece_Hz', 'tied_Hz'}));
figure('Color', 'w');
s.plotMode([], V(:,1), Amplitude = 0.1, Title = 'Mode 1 of the tied halves');

%% 3. Tolerance
% By default Tol is 10% of the median face edge of the finer surface (here
% about 1.5 mm), which absorbs small faceting gaps. Parts that are further
% apart are not tied: a 5 mm gap gives no tied nodes (with a warning)
% unless Tol is raised. Set Tol explicitly when it matters.
fmbd.boxSTL(fullfile(out, 'gap.stl'), [505 0 0], [1000 50 50]);
s2 = fmbd.System();
s2.add(body('left', 25));
s2.add(body('gap', 15));
fprintf('default Tol: %d nodes\n', s2.tie('left', 'gap'));
fprintf('Tol = 6 mm: %d nodes\n', s2.tie('left', 'gap', Tol = 6e-3));

%% 4. Partial bonding with Select
% Select restricts the slave nodes, e.g. to model an adhesive layer that
% covers only part of the interface. Gluing only the lower half (z < 25 mm)
% makes the joint more flexible, mostly for bending in the xz plane.
s2 = fmbd.System();
s2.add(body('left', 25));
s2.add(body('right', 15));
s2.tie('left', 'right', Select = @(X) X(:,3) < 0.025);
s2.fix('left', clamp);
fprintf('half-bonded joint: f1, f2 = %.2f, %.2f Hz (fully bonded %.2f Hz)\n', s2.modes([], 2), f(1));

%% 5. Elastic interfaces
% With Stiffness = [kn kt] (N/m^3, per unit area) the slave nodes are
% connected by springs instead, normal (kn) and tangential (kt) to the
% master surface. This models a compliant layer such as glue, a rubber pad,
% or contact stiffness. A very stiff interface approaches the bonded tie.
kn = 10.^(10:2:16);  f1 = zeros(size(kn));
for k = 1:numel(kn)
    s2 = fmbd.System();
    s2.add(body('left', 25));
    s2.add(body('right', 15));
    s2.tie('left', 'right', Stiffness = kn(k)); % kt = kn
    s2.fix('left', clamp);
    f1(k) = s2.modes([], 1);
end
figure('Color', 'w');
semilogx(kn, f1, 'o-', kn, f(1) + 0*kn, '--');  grid on
xlabel('interface stiffness k_n = k_t (N/m^3)');  ylabel('f_1 (Hz)');
legend('elastic interface', 'bonded tie', 'Location', 'southeast');

%% 6. Ties between elastic and rigid bodies
% A tie between an elastic and a rigid body makes the tied nodes of the
% elastic body a coupling, fixed to the rigid body (whichever order you
% give). Here the right half is rigid: only the left half bends.
s2 = fmbd.System();  s2.add(body('left', 25));
s2.add(fmbd.Body.rigid(fullfile(out, 'right.stl'), 7850, Scale = 1e-3, Name = 'right'));
s2.tie('right', 'left');  s2.fix('left', clamp);
fprintf('rigid right half: f1 = %.2f Hz (elastic %.2f Hz)\n', s2.modes([], 1), f(1));
