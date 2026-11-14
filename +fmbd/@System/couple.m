function name = couple(s, body, sel, opts)
%COUPLE Couple nodes of an elastic body to a new reference point (RBE2).
%   name = s.couple(body, sel) returns the name of the reference point.
%   sel selects nodes by coordinates, @(X) ..., or by index. The first
%   coupling of a body carries its frame.
%   Options: Name, Point (default: centroid of the nodes).
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
