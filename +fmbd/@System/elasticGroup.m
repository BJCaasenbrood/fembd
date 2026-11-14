function B = elasticGroup(s, members)
% Reduced data of the elastic bodies members (tied together): their
% couplings, ties, node supports and node-set input loads, in group
% node numbering.
    go = zeros(1, numel(s.Bodies));  go(members) = s.Uo(members);
    nn = sum(cellfun(@(b) b.nNode, s.Bodies(members)));
    cpl = find(ismember([s.Refs.body], members) & ~cellfun(@isempty, {s.Refs.nodes}));
    C = struct('Name', {}, 'nodes', {}, 'x', {});
    for k = cpl
        C(end+1) = struct('Name', s.Refs(k).Name, 'nodes', go(s.Refs(k).body) + s.Refs(k).nodes(:), ...
                          'x', s.Refs(k).x);                                    %#ok<AGROW>
    end
    Tg = s.Ties(ismember([s.Ties.master], members));
    for t = 1:numel(Tg)
        Tg(t).snodes = go(Tg(t).slave) + Tg(t).snodes;  Tg(t).mnodes = go(Tg(t).master) + Tg(t).mnodes;
    end
    sup = zeros(0, 1);
    for f = s.Fixes(strcmp({s.Fixes.kind}, 'nodes') & ismember([s.Fixes.idx], members))
        d = 3*(go(f.idx) + f.nodes(:)') - 3 + f.dirs(:);  sup = [sup; d(:)];         %#ok<AGROW>
    end
    loads = zeros(3*nn, 0);                   % node-set inputs: static correction shapes
    for g = s.Inputs(strcmp({s.Inputs.kind}, 'body') & ismember([s.Inputs.idx], members))
        f = zeros(3*nn, 1);  f(3*(go(g.idx) + g.nodes(:)) - 3 + g.dir) = 1/numel(g.nodes);
        loads = [loads, f];                                                         %#ok<AGROW>
    end
    name = strjoin(cellfun(@(b) b.Name, s.Bodies(members), 'UniformOutput', false), ' + ');
    B = reduceGroup(s.Bodies(members), C, Tg, unique(sup), loads, ...
                    sum(cellfun(@(b) b.Modes, s.Bodies(members))), name);
    B.cpl = cpl;
end
