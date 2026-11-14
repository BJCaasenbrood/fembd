function name = transmission(s, ids, coeffs, opts)
%TRANSMISSION Belt, gear or cable between coordinates.
%   s.transmission(ids, c) acts on the stretch delta = sum(c_i s_i) of joint
%   coordinates or reference-point coordinates (Dirs 1:6, absolute).
%   With Stiffness it is a spring (plus Damping), without it delta = 0.
%   Example, motor angle N times the link angle:
%     s.transmission({'motor', 'lower'}, [1/N -1], Dirs = {1, 6}, Stiffness = k)
%   Options: Dirs, Stiffness, Damping, Name.
    arguments
        s
        ids (1,:) cell
        coeffs (1,:) double
        opts.Dirs = []
        opts.Stiffness double = []
        opts.Damping (1,1) double = 0
        opts.Name (1,:) char = ''
    end
    nj = numel(ids);
    assert(numel(coeffs) == nj, 'fmbd:transmission', 'One coefficient per entry.');
    dirs = ones(1, nj);
    if iscell(opts.Dirs), dirs = [opts.Dirs{:}];
    elseif ~isempty(opts.Dirs), dirs = opts.Dirs(:)'; end
    assert(numel(dirs) == nj, 'fmbd:transmission', 'One entry of Dirs per coordinate.');
    terms = struct('kind', {}, 'idx', {}, 'dir', {});
    for i = 1:nj
        [kind, k] = lookup(s, ids{i}, true);
        assert(ismember(kind, {'joint', 'ref'}), 'fmbd:transmission', ...
            '%s is not a joint, rigid body or coupling.', string(ids{i}));
        nd = 6;  if strcmp(kind, 'joint'), nd = nnz(s.Joints(k).free); end
        assert(dirs(i) >= 1 && dirs(i) <= nd, 'fmbd:transmission', ...
            '%s has %d coordinate(s).', string(ids{i}), nd);
        terms(i) = struct('kind', kind, 'idx', k, 'dir', dirs(i));
    end
    assert(isempty(opts.Stiffness) || isscalar(opts.Stiffness), 'fmbd:transmission', ...
        'Stiffness must be a scalar (or empty for a rigid transmission).');
    assert(~isempty(opts.Stiffness) || opts.Damping == 0, 'fmbd:transmission', ...
        'Damping needs an elastic transmission (give Stiffness, possibly 0).');
    name = opts.Name;
    if isempty(name), name = sprintf('transmission%d', numel(s.Transmissions) + 1); end
    s.Transmissions(end+1) = struct('Name', name, 'terms', terms, 'coef', coeffs, ...
                                    'k', opts.Stiffness, 'damp', opts.Damping);
    s.built = false;
end
