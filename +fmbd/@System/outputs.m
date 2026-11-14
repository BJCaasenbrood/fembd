function y = outputs(s, q, nu, a)
% Output signals as a struct (fields in the order of OutputName), with
% the accelerations a = dnu/dt for acceleration outputs.
    v = zeros(numel(s.Out), 1);
    for i = 1:numel(s.Out)
        sg = s.Out(i);  [c, Jc] = s.coord(sg.kind, sg.idx, q);
        switch sg.type
            case 'displacement', v(i) = c(sg.dir);
            case 'velocity',     v(i) = Jc(sg.dir, :)*nu;
            otherwise                        % J a + dJ/dt nu
                v(i) = Jc(sg.dir, :)*a;
                if any(nu)
                    e = 1e-6/norm(nu);
                    [~, Jp] = s.coord(sg.kind, sg.idx, s.retract(q, e*nu));
                    [~, Jm] = s.coord(sg.kind, sg.idx, s.retract(q, -e*nu));
                    v(i) = v(i) + (Jp(sg.dir, :) - Jm(sg.dir, :))*nu/(2*e);
                end
        end
    end
    y = cell2struct(num2cell(v), matlab.lang.makeValidName(s.OutputName(:)), 1);
    if isempty(s.Out), y = struct(); end
end
