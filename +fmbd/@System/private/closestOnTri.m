function [d, l] = closestOnTri(p, A, B, C)
% Distance from p to triangles (A, B, C) and the closest-point coordinates.

e1 = B - A;
e2 = C - A;
v = p - A;
a = dot(e1, e1, 2);
b = dot(e1, e2, 2);
c = dot(e2, e2, 2);
r1 = dot(e1, v, 2);  r2 = dot(e2, v, 2);
s = (c.*r1 - b.*r2)./(a.*c - b.^2);
t = (a.*r2 - b.*r1)./(a.*c - b.^2);
l = [1 - s - t, s, t];

out = any(l < 0, 2);
if any(out) % closest point on an edge
    clamp = @(x) min(max(x, 0), 1);
    e3 = C(out,:) - B(out,:);
    u1 = clamp(r1(out)./a(out));
    u2 = clamp(r2(out)./c(out));
    u3 = clamp(dot(p - B(out,:), e3, 2)./dot(e3, e3, 2));
    z = zeros(nnz(out), 1);
    cand = cat(3, [1-u1, u1, z], [1-u2, z, u2], [z, 1-u3, u3]);
    q = @(k) cand(:,1,k).*A(out,:) + cand(:,2,k).*B(out,:) + cand(:,3,k).*C(out,:);
    dk = [vecnorm(p - q(1), 2, 2), vecnorm(p - q(2), 2, 2), vecnorm(p - q(3), 2, 2)];
    [~, k] = min(dk, [], 2);
    m = numel(k);
    l(out,:) = cand((1:m)' + (k - 1)*3*m + [0, m, 2*m]);
end
    d = vecnorm(p - (l(:,1).*A + l(:,2).*B + l(:,3).*C), 2, 2);
end
