%% Tutorial 3 - Materials, body matrices and mass properties
% This tutorial covers:
%
% * isotropic and anisotropic (orthotropic) materials
% * body matrices: stiffness K and consistent mass M of a meshed body
% * mass, centre of mass and inertia (|massProperties|)
% * material damping (|Damping|, Rayleigh)

addpath(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(tempdir, 'fmbdTutorials');  if ~isfolder(out), mkdir(out); end
file = fullfile(out, 'beam.stl');
fmbd.boxSTL(file, [0 0 0], [1000 50 50]);
clamp = @(X) X(:,1) < 1e-6;

%% 1. Isotropic material
% A material is a struct with a 6x6 stiffness D in Voigt notation, order
% [xx yy zz yz xz xy] with engineering shear strains, and a density rho.
steel = fmbd.material(210e9, 0.3, 7850);
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
    b = fmbd.Body.elastic(file, mat{1}, Size = 25, Scale = 1e-3, Modes = 1);
    s = fmbd.System();  s.add(b);  s.fix(b, clamp);
    C = inv(mat{1}.D);  Ex = 1/C(1,1);            % Young's modulus along x
    fEB = 1.875104^2/(2*pi)*sqrt(Ex*0.05^4/12/(mat{1}.rho*0.05^2));
    fprintf('f1 = %6.2f Hz (Euler-Bernoulli %6.2f Hz)\n', s.modes([], 1), fEB);
end

%% 3. Body matrices
% |body.matrices()| returns the body's own K and M (3 dofs per node,
% interleaved [ux1 uy1 uz1 ux2 ...], unconstrained). A free body can move
% rigidly without deforming: K*phi = 0 for its six rigid-body modes.
b = fmbd.Body.elastic(file, steel, Size = 25, Scale = 1e-3);
[Kb, Mb] = b.matrices();
r = b.X - b.Ref;  n = b.nNode;  Phi = zeros(3*n, 6);
for d = 1:3, Phi(d:3:end, d) = 1; end          % translations
Phi(2:3:end, 4) = -r(:,3);  Phi(3:3:end, 4) =  r(:,2);   % rotation about x
Phi(1:3:end, 5) =  r(:,3);  Phi(3:3:end, 5) = -r(:,1);   % rotation about y
Phi(1:3:end, 6) = -r(:,2);  Phi(2:3:end, 6) =  r(:,1);   % rotation about z
fprintf('rigid-body strain energy: %.1e (relative)\n', norm(Phi'*Kb*Phi)/norm(Kb, 1));

%% 4. Mass properties
% |massProperties| gives the mass, centre of mass and inertia (about the
% centre of mass) from the consistent mass matrix; they are exact for the
% mesh. In a system, an elastic body's inertia (and how it changes as the
% body deforms) follows from the same mass matrix.
[m, c, J] = b.massProperties();
mx = 7850*1*0.05*0.05;
fprintf('mass    %.4f kg (exact %.4f), centre %s m\n', m, mx, mat2str(c, 3));
fprintf('inertia %s kg m^2 (exact %s)\n', mat2str(diag(J)', 4), ...
        mat2str(mx*[2*0.05^2, 1 + 0.05^2, 1 + 0.05^2]/12, 4));

%% 5. Material damping
% Damping = [alpha beta] adds Rayleigh damping D = alpha*M + beta*K to an
% elastic body. beta = 2*zeta/omega gives the damping ratio zeta at omega,
% rising linearly with frequency; alpha also damps rigid motion (like a
% drag), so it is usually 0. The modal damping ratio of mode k is then
% beta*omega_k/2, here 1% at the first mode.
w1 = 1.875104^2*sqrt(210e9*0.05^4/12/(7850*0.05^2));  beta = 2*0.01/w1;   % steel, mode 1
b = fmbd.Body.elastic(file, steel, Size = 25, Scale = 1e-3, Modes = 4, Damping = [0 beta]);
s = fmbd.System();  s.add(b);  s.fix(b, clamp);
sys = s.statespace();  p = eig(sys.A);  p = p(imag(p) > 0);  [~, i] = sort(abs(p));
fprintf('damping ratios of the first modes: %s (beta*w/2: %s)\n', ...
        mat2str(-real(p(i(1:2)))'./abs(p(i(1:2)))', 3), mat2str(beta*abs(p(i(1:2)))'/2, 3));
