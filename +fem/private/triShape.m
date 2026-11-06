function N = triShape(l, nf)
%TRISHAPE Tri3/Tri6 shape functions at barycentric points l (m x 3).
%   Tri6 edge nodes sit on edges (1,2) (2,3) (3,1).

if nf == 3
    N = l;
else
    N = [l.*(2*l - 1), 4*l(:,1).*l(:,2), 4*l(:,2).*l(:,3), 4*l(:,3).*l(:,1)];
end
end
