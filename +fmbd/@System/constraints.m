function [g, G] = constraints(s, q, reduce)
% Constraint values g(q) = 0 and their Jacobian G (rows s.keep).
    if nargin < 3, reduce = true; end
    g = zeros(0, 1);  G = zeros(0, s.N);
    for j = 1:numel(s.Jn)
        [c, Jc, gr, Gr] = s.jointKin(j, q);  lt = ~s.Jn(j).free(1:3);
        g = [g; c(lt); gr];  G = [G; Jc(lt, :); Gr];                      %#ok<AGROW>
    end
    for T = s.Tr(arrayfun(@(T) isempty(T.k), s.Tr))
        [d, Jd] = s.stretch(T, q);  g = [g; d];  G = [G; Jd];             %#ok<AGROW>
    end
    for L = s.Lk
        [v, Jv] = s.coord('joint', L.joint, q);
        g = [g; v(L.dir) - L.value];  G = [G; Jv(L.dir, :)];              %#ok<AGROW>
    end
    if reduce, g = g(s.keep);  G = G(s.keep, :); end
end
