function [m, c, J] = massProperties(b)
%MASSPROPERTIES Mass m, centre of mass c (1x3) and inertia J (3x3, about
%   c) of the body: from setMass (or a mesh-free rigid body), otherwise
%   from the consistent mass matrix of the mesh (exact for the mesh).
    if ~isempty(b.MassProps)
        p = b.MassProps;
        m = p.m;
        c = p.c(:)';
        J = p.J;
        return
    end
    [~, M] = b.matrices(MassOnly = true);  Ms = M(1:3:end, 1:3:end);
    m = full(sum(Ms(:)));  c = full(sum(Ms, 1)*b.X)/m;
    Xc = b.X - c;  S = full(Xc'*Ms*Xc);
    J = trace(S)*eye(3) - S;  J = (J + J')/2;
end
