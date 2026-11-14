function name = couple(s, body, sel, opts)
%COUPLE Rigid coupling (RBE2): selected nodes of an elastic body move
%   with a new reference point. Returns the coupling name. The first
%   coupling carries the body's frame. Options: Name, Point.
    arguments
        s
        body
        sel
        opts.Name (1,:) char = ''
        opts.Point double = []
    end
    [kind, i] = s.lookup(body, false);
    assert(strcmp(kind, 'body') && ~s.Bodies{i}.Rigid, 'fmbd:couple', ...
        'Couplings are for elastic bodies; a rigid body is already a reference point.');
    b = s.Bodies{i};  nodes = b.nodes(sel);
    assert(~isempty(nodes), 'fmbd:couple', 'No nodes selected.');
    x = opts.Point;  if isempty(x), x = mean(b.X(nodes, :), 1); end
    name = opts.Name;
    if isempty(name), name = sprintf('%s_ref%d', b.Name, numel(s.Refs) + 1); end
    assert(~ismember(name, s.names()), 'fmbd:couple', 'The name %s is already used.', name);
    s.Refs(end+1) = struct('Name', name, 'body', i, 'nodes', nodes, 'x', x(:)');
    s.built = false;
end
