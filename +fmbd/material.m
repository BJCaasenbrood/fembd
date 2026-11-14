function mat = material(E, nu, rho)
%MATERIAL Isotropic linear-elastic material.
%   mat = fmbd.material(E, nu, rho) returns a struct with E, nu, rho, and D.
%   Voigt order is [xx yy zz yz xz xy] with engineering shear strains.

lam = E*nu/((1 + nu)*(1 - 2*nu));
mu  = E/(2*(1 + nu));
D   = blkdiag(lam*ones(3) + 2*mu*eye(3), mu*eye(3));
mat = struct('E', E, 'nu', nu, 'rho', rho, 'D', D); % assemble the material struct
end
