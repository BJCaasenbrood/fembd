%DEMOPARALLELFLEXURE Linear parallel leaf-spring guide with rigid bodies.
%   A rigid base and a rigid shuttle are connected by two thin elastic
%   leaves. The shuttle moves in x and is stiff in every other direction.
%   The 6x6 stiffness at the shuttle and the first mode are compared with
%   clamped-guided beam theory. Requires Gmsh (set GMSH if needed).
%
%         shuttle (rigid)
%       +================+
%       ||              ||   leaves: L x b x t, elastic
%       ||              ||
%   ====+================+====  base (rigid, fixed)

addpath(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(tempdir, 'femFlexure');  if ~isfolder(out), mkdir(out); end

% geometry (mm)
L = 40;  b = 10;  t = 0.1;  s = 20;            % leaf length, width, thickness; half spacing
fem.boxSTL(fullfile(out, 'base.stl'),    [-s-5 -b/2 -10],   [s+5 b/2 0]);
fem.boxSTL(fullfile(out, 'shuttle.stl'), [-s-5 -b/2 L],     [s+5 b/2 L+10]);
fem.boxSTL(fullfile(out, 'leaf1.stl'),   [-s-t/2 -b/2 0],   [-s+t/2 b/2 L]);
fem.boxSTL(fullfile(out, 'leaf2.stl'),   [ s-t/2 -b/2 0],   [ s+t/2 b/2 L]);

% bodies: rigid parts need only a coarse mesh and a density
steel = fem.material(210e9, 0.3, 7850);
mesh  = {'Scale', 1e-3};
base    = fem.Body(fullfile(out, 'base.stl'),    steel, mesh{:}, Size = 5, Order = 1, Rigid = true);
shuttle = fem.Body(fullfile(out, 'shuttle.stl'), steel, mesh{:}, Size = 5, Order = 1, Rigid = true);
leaf1   = fem.Body(fullfile(out, 'leaf1.stl'),   steel, mesh{:}, Size = 0.8);
leaf2   = fem.Body(fullfile(out, 'leaf2.stl'),   steel, mesh{:}, Size = 0.8);

a = fem.Assembly();
cellfun(@(x) a.add(x), {base, shuttle, leaf1, leaf2});
for leaf = {'leaf1', 'leaf2'}
    a.tie('base', leaf{1});                  % leaves clamped into the rigid parts
    a.tie('shuttle', leaf{1});
end
a.fix('base', [], 1:6);                      % ground the base

[K, M, T] = a.system();
fprintf('%d dofs, %d after constraints\n', a.nDof, size(K, 1));

% 6x6 stiffness at the shuttle centroid: static condensation (Guyan) onto
% its reference dofs [ux uy uz rx ry rz]
ib = a.qdof(a.dof('shuttle'));
Ks = fem.craigBampton(K, M, ib, 0);
disp(array2table(diag(Ks)', 'VariableNames', ...
     {'kx_N_m', 'ky_N_m', 'kz_N_m', 'krx_Nm', 'kry_Nm', 'krz_Nm'}));

% first modes
[V, f] = fem.modes(K, M, 4);

% clamped-guided beam theory (plate value E/(1-nu^2) as upper bound)
I  = b*t^3/12*1e-12;  Lm = L*1e-3;
kx = 2*12*steel.E*I/Lm^3*[1, 1/(1 - steel.nu^2)];
m  = steel.rho*(2*s + 10)*b*10*1e-9 + 2*13/35*steel.rho*L*b*t*1e-9;  % shuttle + leaves
fprintf('kx: FE %.0f N/m, beam theory %.0f - %.0f N/m\n', Ks(1,1), kx);
fprintf('f1: FE %.1f Hz, beam theory %.1f - %.1f Hz\n', f(1), sqrt(kx/m)/(2*pi));
fprintf('support stiffness ratio kz/kx = %.0f, ky/kx = %.0f\n', Ks(3,3)/Ks(1,1), Ks(2,2)/Ks(1,1));

figure('Color', 'w');
for k = 1:4
    subplot(1, 4, k);
    u = T*V(:,k);
    a.plot(u, 0.005/max(abs(u(1:3*size(a.X, 1)))));
    title(sprintf('Mode %d: %.0f Hz', k, f(k)));  axis off
end
