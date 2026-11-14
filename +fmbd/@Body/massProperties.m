function [m, c, J] = massProperties(b)
%MASSPROPERTIES Mass, centre of mass and inertia (about the centre of mass).
%   [m, c, J] = b.massProperties()
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
