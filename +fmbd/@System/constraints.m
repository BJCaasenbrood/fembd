function [g, G] = constraints(s, q, reduce)
% Constraint values g(q) = 0 and their Jacobian G (rows s.keepRows).
    if nargin < 3, reduce = true; end
    g = zeros(0, 1);  G = zeros(0, s.nv);
    for j = 1:numel(s.jointFrames)
        [c, Jc, gr, Gr] = s.jointKin(j, q);  lt = ~s.jointFrames(j).free(1:3);
        g = [g; c(lt); gr];  G = [G; Jc(lt, :); Gr]; %#ok<AGROW>
    end
    for T = s.transData(arrayfun(@(T) isempty(T.k), s.transData))
        [d, Jd] = s.stretch(T, q);
        g = [g; d]; %#ok<AGROW>
        G = [G; Jd]; %#ok<AGROW>
    end
    for L = s.locks
        [v, Jv] = s.coord('joint', L.joint, q);
        g = [g; v(L.dir) - L.value];  G = [G; Jv(L.dir, :)]; %#ok<AGROW>
    end
    if reduce, g = g(s.keepRows);  G = G(s.keepRows, :); end
end
