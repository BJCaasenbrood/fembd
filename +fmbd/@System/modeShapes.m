function [Phi, X0] = modeShapes(s, body)
%MODESHAPES Nodal mode shapes of an elastic body: [Phi, X0] =
%   s.modeShapes(body) with Phi (3n x ne, interleaved [ux1 uy1 uz1 ...])
%   in the body frame and X0 (n x 3) the nodes relative to the frame
%   origin (the first coupling), so the nodes are at r + R*(X0 + Phi*eta)
%   with r, R and eta from simulate (x.r, x.R, x.e).
    s.build();
    [kind, i] = s.lookup(body, false);
    B = s.Bd(s.Ub(i));
    assert(strcmp(kind, 'body') && B.elastic, 'fmbd:modeShapes', 'Not an elastic body.');
    rows = s.Uo(i) + (1:s.Bodies{i}.nNode);  d = 3*rows - [2; 1; 0];
    Phi = B.Phi(d(:), :);  X0 = B.Xl(rows, :);
end
