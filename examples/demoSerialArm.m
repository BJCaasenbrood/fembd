%DEMOSERIALARM Two-link arm with revolute joints (couplings and joints).
%   Two aluminium links are hinged to ground and to each other about z, with
%   torsion springs in the joints (servo stiffness). The links are elastic,
%   attached to the joints through rigid couplings (RBE2) on their end
%   faces; set rigidLinks = true to make them rigid instead. The arm is
%   linearised at several elbow angles, held by the joint torques, to show
%   how the dynamics depend on the configuration.
%
%   ground --(shoulder)-- link1 --(elbow)-- link2

rigidLinks = false;
addpath(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(tempdir, 'fmbdArm');  if ~isfolder(out), mkdir(out); end

% geometry (mm): links along x, 200 x 20 x 10
fmbd.boxSTL(fullfile(out, 'link1.stl'), [0 -10 -5],   [200 10 5]);
fmbd.boxSTL(fullfile(out, 'link2.stl'), [200 -10 -5], [400 10 5]);
alu = fmbd.material(1e9, 0.33, 2700);
kj = 200; % joint stiffness, Nm/rad

arm = @(elbow) twoLink(out, alu, rigidLinks, kj, elbow);
s = arm(pi/3);
% the joint springs are relaxed at zero, so the servo torques hold the
% elbow at pi/3: linearise about that equilibrium
x0 = s.equilibrium(Input = kj*[0; pi/3]);
[f, V] = s.modes(x0, 6);
fprintf('%d model dofs\n', s.nDof);
disp(table((1:6)', f, 'VariableNames', {'mode', 'f_Hz'}));
figure('Color', 'w');

for k = 1:4
    subplot(2, 2, k);  s.plotMode(x0, V(:,k), Amplitude = 0.03);
    title(sprintf('Mode %d: %.1f Hz', k, f(k)));
    view(2); axis off
end

% actuator-level model: joint torques to joint angles, with 4 modes
s.output('q1', 'shoulder');   s.output('q2', 'elbow'); % (inputs tau1, tau2: see twoLink)
S = s.statespace(x0, Modes = 4);
disp('Joint-space static stiffness (Nm/rad):');  disp(inv(-S.C*(S.A\S.B)));

% configuration dependence: at each elbow angle, the servo torques hold the
% arm against the joint springs (relaxed at zero), and the model is
% linearised about that equilibrium
angles = 0:15:90;  fa = zeros(numel(angles), 2);
for k = 1:numel(angles)
    s = arm(deg2rad(angles(k)));
    x = s.equilibrium(Input = kj*[0; deg2rad(angles(k))]);
    fk = s.modes(x, 2);  fa(k,:) = fk(1:2);
end
disp(table(angles', fa(:,1), fa(:,2), 'VariableNames', {'elbow_deg', 'mode1_Hz', 'mode2_Hz'}));
s = arm(deg2rad(60));  x = s.equilibrium(Input = kj*[0; deg2rad(60)]);
[f1, V1] = s.modes(x, 1);
figure('Color', 'w');  s.plotMode(x, V1, Amplitude = 0.03);
title(sprintf('Elbow 60^o: mode 1 at %.1f Hz', f1));  view(2);

function s = twoLink(out, mat, rigid, kj, theta)
% Two-link arm with the elbow at a given angle (Position).
s = fmbd.System();
opt = {'Scale', 1e-3, 'Size', 5, 'Rigid', rigid, 'Modes', 6};
s.add(fmbd.Body(fullfile(out, 'link1.stl'), mat, opt{:}, Name = 'link1'));
s.add(fmbd.Body(fullfile(out, 'link2.stl'), mat, opt{:}, Name = 'link2'));
if rigid % joint attachment: the body, or an end-face coupling
    [p1, p2a, p2b] = deal('link1', 'link1', 'link2');
else
    p1  = s.couple('link1', @(X) X(:,1) < 1e-6); % hub of link 1
    p2a = s.couple('link1', @(X) X(:,1) > 0.2 - 1e-6); % elbow, link 1 side
    p2b = s.couple('link2', @(X) X(:,1) < 0.2 + 1e-6); % elbow, link 2 side
end
s.joint('ground', p1, 'revolute', Point = [0 0 0],   Stiffness = kj, Name = 'shoulder');
s.joint(p2a, p2b,     'revolute', Point = [0.2 0 0], Stiffness = kj, Name = 'elbow', ...
        Position = theta);
s.input('tau1', 'shoulder'); 
s.input('tau2', 'elbow');
end
