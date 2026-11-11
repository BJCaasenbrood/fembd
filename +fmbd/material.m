function mat = material(E, nu, rho)
%MATERIAL Isotropic linear-elastic material.
%   mat = fmbd.material(E, nu, rho) returns a struct with Young's modulus E,
%   Poisson ratio nu, density rho and the 6x6 Voigt stiffness D. Voigt order
%   is [xx yy zz yz xz xy] with engineering shear strains.
%
%   Any other linear material can be used by building the struct directly,
%   e.g. mat = struct('D', D, 'rho', rho) for an anisotropic D.

lam = E*nu/((1 + nu)*(1 - 2*nu));
mu  = E/(2*(1 + nu));
D   = blkdiag(lam*ones(3) + 2*mu*eye(3), mu*eye(3));
mat = struct('E', E, 'nu', nu, 'rho', rho, 'D', D);
end
