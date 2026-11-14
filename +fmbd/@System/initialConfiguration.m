function q = initialConfiguration(s, q)
% Satisfy the constraints and the joint positions exactly (hard), with
% elastic transmissions unstretched unless all their coordinates are
% set; the other joint coordinates stay as close to zero as possible
% (soft, in the null space of the hard equations). Positions are ramped
% up for large angles; elastic bodies stay undeformed.
    tgt = zeros(0, 3);  soft = zeros(0, 2);  % [joint, coordinate, value], [joint, coordinate]
    for j = 1:numel(s.Joints)
        v = s.Joints(j).position;  nf = nnz(s.Joints(j).free);
        for d = 1:numel(v), tgt(end+1, :) = [j, d, v(d)]; end                 %#ok<AGROW>
        if isempty(v), soft = [soft; j*ones(nf, 1), (1:nf)']; end              %#ok<AGROW>
    end
    given = @(T) all(arrayfun(@(t) strcmp(t.kind, 'joint') && ...
                 ~isempty(s.Joints(t.idx).position), T.terms));
    relax = s.Tr(arrayfun(@(T) ~isempty(T.k) && ~given(T), s.Tr));
    rigid = [s.Bd.cols];  ec = [s.Bd.ecols];
    if ~isempty(ec)                          % keep the modal coordinates at zero
        el = arrayfun(@(B) B.cols(7:end), s.Bd, 'UniformOutput', false);
        rigid = setdiff(rigid, [el{:}]);
    end
    ns = max(1, ceil(max([abs(tgt(:,3)); 0])/0.3));
    for step = 1:ns
        for it = 1:50
            [g, G] = s.constraints(q, false);
            for i = 1:size(tgt, 1)
                [v, Jv] = s.coord('joint', tgt(i,1), q);
                g = [g; v(tgt(i,2)) - step/ns*tgt(i,3)];  G = [G; Jv(tgt(i,2), :)];   %#ok<AGROW>
            end
            for T = relax
                [d, Jd] = s.stretch(T, q);  g = [g; d];  G = [G; Jd];               %#ok<AGROW>
            end
            gs = zeros(size(soft, 1), 1);  Gs = zeros(size(soft, 1), s.N);
            for i = 1:size(soft, 1)
                [v, Jv] = s.coord('joint', soft(i,1), q);  gs(i) = v(soft(i,2));  Gs(i,:) = Jv(soft(i,2), :);
            end
            Gh = G(:, rigid);  dh = -lsqminnorm(Gh, g);   % hard equations, minimum norm
            if ~isempty(gs)                               % soft: in the null space of the hard ones
                Nh = null(full(Gh));  if isempty(Gh), Nh = eye(numel(rigid)); end
                dh = dh - Nh*lsqminnorm(Gs(:, rigid)*Nh, gs + Gs(:, rigid)*dh);
            end
            if max(abs([g; 0])) < 1e-12 && norm(dh) < 1e-12, break, end
            dq = zeros(s.N, 1);  dq(rigid) = dh;
            q = s.retract(q, dq);
        end
    end
    assert(isempty(g) || max(abs(g)) < 1e-9, 'fmbd:System', ['The initial configuration ' ...
        'cannot satisfy the joints, transmissions and joint positions (residual %.1e).'], max(abs(g)));
end
