%DEMOPARALLELFLEXURE Linear parallel leaf-spring guide with rigid bodies.
%   A rigid base and a rigid shuttle are connected by two thin elastic
%   leaves. The shuttle moves in x and is stiff in every other direction.
%   The 6x6 stiffness at the shuttle and the first mode are compared with
%   clamped-guided beam theory. Requires Gmsh (set GMSH if needed).
%
%         shuttle (rigid)
%       +================+
%       ||              ||  leafsprings: L-b-t, (elastic)
%       ||              ||
%   ====+================+====  base (rigid)

addpath(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(tempdir, 'fmbdFlexure');  if ~isfolder(out), mkdir(out); end

% geometry (mm)
L = 40;  b = 10;  t = 0.1;  h = 20; % leaf length, width, thickness; half spacing
fmbd.boxSTL(fullfile(out, 'base.stl'),    [-h-5 -b/2 -10],   [h+5 b/2 0]);
fmbd.boxSTL(fullfile(out, 'shuttle.stl'), [-h-5 -b/2 L],     [h+5 b/2 L+10]);
fmbd.boxSTL(fullfile(out, 'leaf1.stl'),   [-h-t/2 -b/2 0],   [-h+t/2 b/2 L]);
fmbd.boxSTL(fullfile(out, 'leaf2.stl'),   [ h-t/2 -b/2 0],   [ h+t/2 b/2 L]);

% bodies: the rigid parts from their STL surfaces, the leaves meshed
steel = fmbd.material(210e9, 0.3, 7850);
s = fmbd.System();
s.add(fmbd.Body.rigid(fullfile(out, 'base.stl'),    steel.rho, Scale = 1e-3, Name = 'base'));
s.add(fmbd.Body.rigid(fullfile(out, 'shuttle.stl'), steel.rho, Scale = 1e-3, Name = 'shuttle'));
for leaf = {'leaf1', 'leaf2'}
    s.add(fmbd.Body.elastic(fullfile(out, [leaf{1} '.stl']), steel, Scale = 1e-3, Size = 0.8, ...
                            Modes = 6, Name = leaf{1}, Order=1));
    s.tie('base', leaf{1});                    % leaves clamped into the rigid parts
    s.tie('shuttle', leaf{1});
end
s.fix('base');                                 % ground the base

% 6x6 stiffness at the shuttle: the static gain from forces and moments on
% the shuttle to its displacements and rotations (kept by the static
% correction), inverted
s.input('F', 'shuttle');  s.output('u', 'shuttle');
S = s.statespace();
Ks = inv(-S.C*(S.A\S.B));
disp(array2table(diag(Ks)', 'VariableNames', ...
     {'kx_N_m', 'ky_N_m', 'kz_N_m', 'krx_Nm', 'kry_Nm', 'krz_Nm'}));

% first modes vs clamped-guided beam theory (plate value E/(1-nu^2) as upper bound)
[f, V] = s.modes([], 4);
I  = b*t^3/12*1e-12;  Lm = L*1e-3;
kx = 2*12*steel.E*I/Lm^3*[1, 1/(1 - steel.nu^2)];
m  = steel.rho*(2*h + 10)*b*10*1e-9 + 2*13/35*steel.rho*L*b*t*1e-9;  % shuttle + leaves
fprintf('kx: %.0f N/m, beam theory %.0f - %.0f N/m\n', Ks(1,1), kx);
fprintf('f1: %.1f Hz, beam theory %.1f - %.1f Hz\n', f(1), sqrt(kx/m)/(2*pi));
fprintf('support stiffness ratio kz/kx = %.0f, ky/kx = %.0f\n', Ks(3,3)/Ks(1,1), Ks(2,2)/Ks(1,1));
figure('Color', 'w');
for k = 1:4
    subplot(1, 4, k);  s.plotMode([], V(:,k), Amplitude = 0.005);
    title(sprintf('Mode %d: %.0f Hz', k, f(k)));  axis off
end
