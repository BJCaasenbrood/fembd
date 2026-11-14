function name = add(s, b)
%ADD Add a body to the system.
%   s.add(b) adds fmbd.Body b and returns its name. A rigid body is also a
%   reference point with the same name.
    arguments
        s
        b (1,1) fmbd.Body
    end
    assert(~any(cellfun(@(c) c == b, s.Bodies)), 'fmbd:add', '%s is already added.', b.Name);
    assert(~ismember(b.Name, s.names()), 'fmbd:add', 'The name %s is already used.', b.Name);
    s.Bodies{end+1} = b;  name = b.Name;
    if b.Rigid
        s.Refs(end+1) = struct('Name', b.Name, 'body', numel(s.Bodies), 'nodes', [], 'x', b.Ref);
    end
    s.built = false;
end
