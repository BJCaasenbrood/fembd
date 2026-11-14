function sig = signalSpec(s, name, id, sel, dirs, type)
% One signal per direction: name, or name_dir if there are several.
    used = [s.InputName, s.OutputName];
    assert(~any(strcmp(used, name) | startsWith(used, [name '_'])), 'fmbd:io', ...
        'Signal name %s is already used.', name);
    [kind, k] = s.lookup(id, isempty(sel));
    nodes = [];
    if strcmp(kind, 'body')
        nodes = s.Bodies{k}.nodes(sel);  if isempty(sel), nodes = (1:s.Bodies{k}.nNode)'; end
        assert(~isempty(nodes), 'fmbd:io', '%s: no nodes selected.', name);
        if isempty(dirs), dirs = 1:3; end
    elseif isempty(dirs)
        dirs = 1:6;  if strcmp(kind, 'joint'), dirs = 1:nnz(s.Joints(k).free); end
    end
    sig = struct('Name', {}, 'kind', {}, 'idx', {}, 'dir', {}, 'nodes', {}, 'type', {});
    for d = dirs(:)'
        nm = name;  if numel(dirs) > 1, nm = sprintf('%s_%d', name, d); end
        sig(end+1) = struct('Name', nm, 'kind', kind, 'idx', k, 'dir', d, ...
                            'nodes', nodes, 'type', type); %#ok<AGROW>
    end
end
