%% Tutorial 3 - Materials, mass matrices and mass properties
% This tutorial covers:
%
% * isotropic and anisotropic (orthotropic) materials
% * body matrices, and mass and inertia from the mass matrix
% * consistent vs lumped mass

addpath(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(tempdir, 'femTutorials');  if ~isfolder(out), mkdir(out); end
file = fullfile(out, 'beam.stl');
fem.boxSTL(file, [0 0 0], [1000 50 50]);
clampX0 = @(X) X(:,1) < 1e-9;

%% 1. Isotropic material
% A material is a struct with a 6x6 stiffness D in Voigt notation, order
% [xx yy zz yz xz xy] with engineering shear strains, and a density rho.
steel = fem.material(210e9, 0.3, 7850);
disp('D of steel (GPa):');  disp(steel.D/1e9);

%% 2. Orthotropic material
% Any symmetric positive-definite D can be used. For an orthotropic
% material (e.g. a unidirectional composite), build the compliance S from
% the engineering constants and invert it. Here the fibres run along x.
E   = [140e9 10e9 10e9];          % E1 E2 E3
nu  = [0.3 0.3 0.4];              % nu12 nu13 nu23
G   = [3.5e9 5e9 5e9];            % G23 G13 G12 (Voigt order yz xz xy)
S = [ 1/E(1)       -nu(1)/E(1)  -nu(2)/E(1);
     -nu(1)/E(1)    1/E(2)      -nu(3)/E(2);
     -nu(2)/E(1)   -nu(3)/E(2)   1/E(3)];
cfrp = struct('D', inv(blkdiag(S, diag(1./G))), 'rho', 1600);

% Compare a steel and a CFRP cantilever. With the fibres along the beam,
% bending is governed by E1. E1/rho is much higher than for steel, so the
% CFRP beam is faster; its low shear moduli make it slightly softer than
% beam theory predicts.
for mat = {steel, cfrp}
    b = fem.Body(file, mat{1}, Size = 25, Scale = 1e-3);
    a = fem.Assembly();  a.add(b);  a.fix(1, clampX0);
    [K, M] = a.system();  [~, f] = fem.modes(K, M, 1);
    C = inv(mat{1}.D);  Ex = 1/C(1,1);            % Young's modulus along x
    fEB = 1.875104^2/(2*pi)*sqrt(Ex*0.05^4/12/(mat{1}.rho*0.05^2));
    fprintf('f1 = %6.2f Hz (Euler-Bernoulli %6.2f Hz)\n', f, fEB);
end

%% 3. Body matrices
% |body.matrices()| returns the body's own K and M (3 dofs per node,
% unconstrained). A body that is not supported can move rigidly without
% deforming: K*phi = 0 for the six rigid-body modes, and phi'*M*phi gives
% its mass and inertia.
b = fem.Body(file, steel, Size = 25, Scale = 1e-3);
[Kb, Mb] = b.matrices();
c = b.Ref;                                     % centroid
r = b.X - c;  n = b.nNode;  Phi = zeros(3*n, 6);
for d = 1:3, Phi(d:3:end, d) = 1; end          % translations
Phi(2:3:end, 4) = -r(:,3);  Phi(3:3:end, 4) =  r(:,2);   % rotation about x
Phi(1:3:end, 5) =  r(:,3);  Phi(3:3:end, 5) = -r(:,1);   % rotation about y
Phi(1:3:end, 6) = -r(:,2);  Phi(2:3:end, 6) =  r(:,1);   % rotation about z

m = 7850*1*0.05*0.05;
Iref = m*[2*0.05^2, 1 + 0.05^2, 1 + 0.05^2]/12;
fprintf('rigid-body strain energy: %.1e (relative)\n', norm(Phi'*Kb*Phi)/norm(Kb, 1));
fprintf('mass    %.4f kg (exact %.4f)\n', Phi(:,1)'*Mb*Phi(:,1), m);
fprintf('inertia %s kg m^2 (exact %s)\n', mat2str(diag(Phi(:,4:6)'*Mb*Phi(:,4:6))', 4), mat2str(Iref, 4));

%% 4. Consistent vs lumped mass
% The default mass matrix is consistent. With Lumped = true it is
% diagonal, using HRZ lumping: each element's mass diagonal is scaled to
% keep the element mass. This is useful for explicit time integration.
% Note: rigid ties spread slave masses onto master nodes, so the reduced M
% is then no longer diagonal.
[~, Ml] = b.matrices(Lumped = true);
ml = full(diag(Ml));
fprintf('lumped M diagonal: %d, min entry %.2e kg, total %.4f kg\n', isdiag(Ml), min(ml), sum(ml)/3);

a = fem.Assembly();  a.add(b);  a.fix(1, clampX0);
[K, M]  = a.system();               [~, fc] = fem.modes(K, M, 6);
[K, Ml] = a.system(Lumped = true);  [~, fl] = fem.modes(K, Ml, 6);
disp(table(fc, fl, 100*(fl./fc - 1), 'VariableNames', {'consistent_Hz', 'lumped_Hz', 'diff_pct'}));
