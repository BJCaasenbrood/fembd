function sig = compileSignals(s, list)
% Signals -> (kind, idx, dir, type); node sets become average-displacement signals.
    sig = struct('kind', {list.kind}, 'idx', {list.idx}, 'dir', {list.dir}, 'type', {list.type});
    for c = find(strcmp({list.kind}, 'body'))
        i = list(c).idx;
        u = s.bodyUnit(i);
        B = s.units(u);
        nodes = s.bodyOffset(i) + list(c).nodes(:);
        d = 3*nodes' - [2; 1; 0];  P = zeros(3, B.ne);
        if B.ne > 0, P = reshape(mean(reshape(B.Phi(d(:), :), 3, numel(nodes), B.ne), 2), 3, B.ne); end
        s.nodeSets(end+1) = struct('body', u, 'u0', mean(B.Xl(nodes, :), 1)', 'Phi', P, ...
                             'x0', mean(s.Bodies{i}.X(list(c).nodes, :), 1)');
        sig(c).kind = 'node';  sig(c).idx = numel(s.nodeSets);
    end
end
