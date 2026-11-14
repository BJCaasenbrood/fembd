function V = vertices(s, i, x, k, scale)
% Node positions of body i at step k (deformation magnified by scale);
% the centre of mass for a point mass.
    u = s.Ub(i);  B = s.Bd(u);  rows = s.Uo(i) + (1:s.Bodies{i}.nNode);
    if isempty(rows), V = x.r(:,u,k)';  return, end
    U = B.Xl(rows, :);
    if B.elastic
        d = 3*rows - [2; 1; 0];
        U = U + scale*reshape(B.Phi(d(:), :)*x.e(B.ecols, k), 3, [])';
    end
    V = U*x.R(:,:,u,k)' + x.r(:,u,k)';
end
