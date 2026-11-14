function name = transmission(s, ids, coeffs, opts)
%TRANSMISSION Belt, gear or cable ratio between coordinates: the stretch
%   delta = sum_i c_i s_i of joint coordinates or reference-point
%   coordinates (Dirs 1:6 = [ux uy uz rx ry rz], absolute).
%     Stiffness = k   elastic: a spring k on delta (with Damping)
%     (none)          rigid: delta = 0
%   Dirs gives the coordinate of each entry (default 1). Examples (ratio
%   N = theta_motor/theta_driven):
%     s.transmission({'motor', 'shoulder'}, [1/N -1], Stiffness = k)
%     s.transmission({'motor', 'lower'}, [1/N -1], Dirs = {1, 6})   % absolute rz
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
        [kind, k] = s.lookup(ids{i}, true);
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
