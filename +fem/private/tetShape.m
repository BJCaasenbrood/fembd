function [N, dN] = tetShape(xi, nn)
%TETSHAPE Tet4/Tet10 shape functions at reference points xi (ng x 3).
%   N is ng x nn, dN = dN/dxi is nn x 3 x ng. Node ordering follows Gmsh;
%   Tet10 edge nodes sit on edges (1,2) (2,3) (3,1) (4,1) (4,3) (4,2).

ng = size(xi, 1);
L  = [1 - sum(xi, 2), xi];
dL = [-1 -1 -1; eye(3)];
if nn == 4
    N  = L;
    dN = repmat(dL, 1, 1, ng);
else
    E  = [1 2; 2 3; 3 1; 4 1; 4 3; 4 2];
    N  = [L.*(2*L - 1), 4*L(:, E(:,1)).*L(:, E(:,2))];
    Lp = permute(L, [2 3 1]);
    dN = [(4*Lp - 1).*dL;
          4*(Lp(E(:,2),:,:).*dL(E(:,1),:) + Lp(E(:,1),:,:).*dL(E(:,2),:))];
end
end
