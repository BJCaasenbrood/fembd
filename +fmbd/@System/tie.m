function ns = tie(s, master, slave, opts)
%TIE Bond the surface of slave to the surface of master.
%   n = s.tie(master, slave) ties the slave surface nodes within Tol of the
%   master surface and returns their number. Use the finer mesh as slave.
%   Tied elastic bodies are reduced together. A tie to a rigid body makes
%   the tied nodes a coupling with a fixed joint.
%   Options: Tol, Stiffness = [kn kt] (elastic interface, per unit area),
%   Select (selector on the slave nodes).
    arguments
        s
        master
        slave
        opts.Tol double = []
        opts.Stiffness double = []
        opts.Select = []
    end
    [~, im] = s.lookup(master, false);  [~, is] = s.lookup(slave, false);
    bm = s.Bodies{im};  bs = s.Bodies{is};
    assert(~isempty(bm.F) && ~isempty(bs.F), 'fmbd:tie', 'A point mass has no surface; use a joint.');
    assert(~(bm.Rigid && bs.Rigid), 'fmbd:tie', 'Two rigid bodies cannot be tied; use a fixed joint.');
    if bs.Rigid, [bm, bs, im, is] = deal(bs, bm, is, im); end % the rigid body is the master
    assert(~bm.Rigid || isempty(opts.Stiffness), 'fmbd:tie', ['An elastic interface to a ' ...
        'rigid body is not supported; use a coupling and a joint with Stiffness.']);
    [sn, mn, Nw, nrm, area] = tieNodes(bm, bs, opts.Tol, opts.Select);
    ns = numel(sn);
    if ns == 0, warning('fmbd:tie', 'No slave nodes within the tolerance.'); return, end
    if bm.Rigid
        c = s.couple(bs.Name, sn, Name = sprintf('%s_tie%d', bs.Name, numel(s.Refs) + 1));
        s.joint(bm.Name, c, 'fixed', Point = s.Refs(end).x, Name = sprintf('%s_tie', c));
    else
        s.Ties(end+1) = struct('master', im, 'slave', is, 'snodes', sn, 'mnodes', mn, ...
                               'N', Nw, 'n', nrm, 'A', area, 'k', opts.Stiffness);
        s.built = false;
    end
end
