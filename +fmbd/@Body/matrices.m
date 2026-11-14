function [K, M] = matrices(b, opts)
%MATRICES Body stiffness K and mass M (sparse, nDof x nDof); K = 0 for
%   a rigid body (and with MassOnly = true).
%   M is consistent, or diagonal with Lumped = true (HRZ: element mass
%   diagonals scaled to preserve the element mass; row-sum lumping
%   would give negative Tet10 corner masses).
    arguments
        b
        opts.Lumped (1,1) logical = false
        opts.Chunk (1,1) double = 2e4 % elements per vectorised block
        opts.MassOnly (1,1) logical = false
    end
    chunk = opts.Chunk;
    [ne, nn] = size(b.T);
    nd = 3*nn;
    p = 1 + (nn == 10);
    [xk, wk] = tetQuad(2*p - 2);  [~, dNk] = tetShape(xk, nn);
    [xm, wm] = tetQuad(2*p);      [Nm, dNm] = tetShape(xm, nn);
    D = [];  if ~b.Rigid, D = b.Mat.D; end
    K = sparse(b.nDof, b.nDof);  Ms = sparse(b.nNode, b.nNode);
    bad = false(1, ne); % elements with det(J) <= 0
    for e0 = 1:chunk:ne
        Te = b.T(e0:min(e0 + chunk - 1, ne), :)';  nc = size(Te, 2);
        Xe = permute(reshape(b.X(Te, :), nn, nc, 3), [1 3 2]);
        Ke = zeros(nd, nd, nc);  Me = zeros(nn, nn, nc);
        for g = 1:numel(wk)*~(b.Rigid || opts.MassOnly)
            [G, dJ] = gradients(dNk(:,:,g), Xe);
            B = zeros(6, nd, nc);
            B(1, 1:3:end, :) = G(1,:,:);  B(2, 2:3:end, :) = G(2,:,:);
            B(3, 3:3:end, :) = G(3,:,:);
            B(4, 2:3:end, :) = G(3,:,:);  B(4, 3:3:end, :) = G(2,:,:);
            B(5, 1:3:end, :) = G(3,:,:);  B(5, 3:3:end, :) = G(1,:,:);
            B(6, 1:3:end, :) = G(2,:,:);  B(6, 2:3:end, :) = G(1,:,:);
            Ke = Ke + pagemtimes(B, 'transpose', pagemtimes(D, B), 'none').*(wk(g)*dJ);
        end
        for g = 1:numel(wm)
            [~, dJ] = gradients(dNm(:,:,g), Xe);
            bad(e0:e0+nc-1) = bad(e0:e0+nc-1) | dJ(:)' <= 0;
            Me = Me + (Nm(g,:)'*Nm(g,:)).*(wm(g)*dJ);
        end
        dof = reshape(3*permute(Te, [3 1 2]) - [2; 1; 0], nd, nc);
        K  = K  + sparse(repmat(dof, nd, 1), repelem(dof, nd, 1), Ke(:), b.nDof, b.nDof);
        if opts.Lumped
            Me = reshape(Me, nn*nn, nc);  d = Me(1:nn+1:end, :);
            Ms = Ms + sparse(Te, Te, d.*sum(Me, 1)./sum(d, 1), b.nNode, b.nNode);
        else
            Ms = Ms + sparse(repmat(Te, nn, 1), repelem(Te, nn, 1), Me(:), b.nNode, b.nNode);
        end
    end
    if any(bad)
        warning('fmbd:Body', '%s: %d element(s) with det(J) <= 0; results are unreliable.', ...
                b.Name, nnz(bad));
    end
    K = (K + K')/2;
    M = kron((Ms + Ms')/2*b.Mat.rho, speye(3));
    if b.Rigid && ~isempty(b.MassProps) % mass from setMass, not from the mesh
        M = sparse(b.nDof, b.nDof);
    end
end
