function [d, Jd] = stretch(s, T, q)
% Transmission stretch delta = sum c_i s_i.
    d = 0;  Jd = zeros(1, s.N);
    for i = 1:numel(T.terms)
        [v, Jv] = s.coord(T.terms(i).kind, T.terms(i).idx, q);
        d = d + T.coef(i)*v(T.terms(i).dir);  Jd = Jd + T.coef(i)*Jv(T.terms(i).dir, :);
    end
end
