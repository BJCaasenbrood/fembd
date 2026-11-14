function [kind, k] = lookup(s, id, whole)
% Resolve id to a 'joint', 'ref' or 'body' index. Joints and reference
% points take precedence when no node selection is given (whole).
    if isa(id, 'fmbd.Body')
        k = find(cellfun(@(b) b == id, s.Bodies));
        assert(isscalar(k), 'fmbd:System', '%s is not part of the system.', id.Name);
        kind = 'body';
        if whole && id.Rigid, k = find([s.Refs.body] == k & cellfun(@isempty, {s.Refs.nodes}));  kind = 'ref'; end
        return
    end
    if whole
        k = find(strcmp({s.Joints.Name}, id));  kind = 'joint';
        if isscalar(k), return, end
        k = find(strcmp({s.Refs.Name}, id));    kind = 'ref';
        if isscalar(k), return, end
    end
    k = find(strcmp(cellfun(@(b) b.Name, s.Bodies, 'UniformOutput', false), id));
    kind = 'body';
    assert(isscalar(k), 'fmbd:System', 'Unknown name %s.', string(id));
    if whole && s.Bodies{k}.Rigid, [kind, k] = s.lookup(s.Bodies{k}, true); end
end
