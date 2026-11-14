function addJoint(s, iA, iB, p, R, free, k, c)
% Joint between reference points iA (0 = ground) and iB at point p with
% frame R (world, defining configuration; anchors are world-aligned).
    aA = 0;  a = p;
    if iA > 0, aA = s.refFrames(iA).anchor;  a = p - s.refFrames(iA).x0 + s.refFrames(iA).s; end
    aB = s.refFrames(iB).anchor;  b = p - s.refFrames(iB).x0 + s.refFrames(iB).s;
    assert(aA ~= aB, 'fmbd:joint', 'A joint connects a reference point to itself.');
    kk = zeros(6, 1);  cc = zeros(6, 1);
    kk(free) = k(:).*ones(nnz(free), 1);  cc(free) = c(:).*ones(nnz(free), 1);
    s.jointFrames(end+1) = struct('A', aA, 'B', aB, 'a', a, 'b', b, 'QA', R, 'QB', R, ...
                         'free', logical(free), 'k', kk, 'c', cc, 'pairs', rotationPairs(free));
end
