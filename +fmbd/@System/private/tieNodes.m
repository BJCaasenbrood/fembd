function [sn, mn, N, n, A] = tieNodes(bm, bs, tol, select)
% Slave surface nodes sn within tol of the master surface, their master
% face nodes mn, shape-function weights N, master normals n and slave
% nodal areas A.
Fm = bm.F;
Am = bm.X(Fm(:,1),:);
Bm = bm.X(Fm(:,2),:);
Cm = bm.X(Fm(:,3),:);
if isempty(tol)
    edge = @(b, F) median(vecnorm(b.X(F(:,[2 3 1]),:) - b.X(F(:,1:3),:), 2, 2));
    tol = 0.1*min(edge(bm, Fm), edge(bs, bs.F));
end
sn = unique(bs.F(:));
if ~isempty(select), sn = sn(select(bs.X(sn,:))); end
lo = min(min(Am, Bm), Cm) - tol;  hi = max(max(Am, Bm), Cm) + tol;
sn = sn(all(bs.X(sn,:) >= min(lo) & bs.X(sn,:) <= max(hi), 2));
d = inf(numel(sn), 1);
f = zeros(numel(sn), 1);
l = zeros(numel(sn), 3);
for i = 1:numel(sn)
    p = bs.X(sn(i),:);
    c = find(all(p >= lo & p <= hi, 2));
    if isempty(c), continue, end
    [dc, lc] = closestOnTri(p, Am(c,:), Bm(c,:), Cm(c,:));
    [d(i), k] = min(dc);
    f(i) = c(k);
    l(i,:) = lc(k,:);
end
on = d <= tol;
sn = sn(on);
f = f(on);
l = l(on,:);
mn = Fm(f, :);  N = triShape(l, size(Fm, 2));
n = cross(Bm(f,:) - Am(f,:), Cm(f,:) - Am(f,:), 2);  n = n./vecnorm(n, 2, 2);
% nodal areas of the slave interface (face weights) from slave faces that
% touch a tied node and face the master surface
Fs = bs.F;
nf = cross(bs.X(Fs(:,2),:) - bs.X(Fs(:,1),:), bs.X(Fs(:,3),:) - bs.X(Fs(:,1),:), 2);
nm = zeros(bs.nNode, 3);  nm(sn,:) = n;
nm = reshape(sum(reshape(nm(Fs,:), [size(Fs) 3]), 2), [], 3);
on = dot(nf, nm, 2) < -0.5*vecnorm(nf, 2, 2).*vecnorm(nm, 2, 2);
Fs = Fs(on,:);  ar = vecnorm(nf(on,:), 2, 2)/2;
w = [1 1 1]/3;  if size(Fs, 2) == 6, w = [3 3 3 16 16 16]/57; end
nodeA = accumarray(Fs(:), reshape(ar.*w, [], 1), [bs.nNode 1]);
A = nodeA(sn);
end
