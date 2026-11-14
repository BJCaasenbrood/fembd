function y = outputs(s, q, nu, a)
% Output signals as a struct (fields in the order of OutputName), with
% the accelerations a = dnu/dt for acceleration outputs.
    v = zeros(numel(s.outSignals), 1);
    for i = 1:numel(s.outSignals)
        sg = s.outSignals(i);  [c, Jc] = coord(s, sg.kind, sg.idx, q);
        switch sg.type
            case 'displacement', v(i) = c(sg.dir);
            case 'velocity',     v(i) = Jc(sg.dir, :)*nu;
            otherwise % J a + dJ/dt nu
                v(i) = Jc(sg.dir, :)*a;
                if any(nu)
                    e = 1e-6/norm(nu);
                    [~, Jp] = coord(s, sg.kind, sg.idx, retract(s, q, e*nu));
                    [~, Jm] = coord(s, sg.kind, sg.idx, retract(s, q, -e*nu));
                    v(i) = v(i) + (Jp(sg.dir, :) - Jm(sg.dir, :))*nu/(2*e);
                end
        end
    end
    y = cell2struct(num2cell(v), matlab.lang.makeValidName(s.OutputName(:)), 1);
    if isempty(s.outSignals), y = struct(); end
end
