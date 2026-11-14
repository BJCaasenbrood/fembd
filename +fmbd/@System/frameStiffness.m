function K = frameStiffness(s, q, u, lam)
% Tangent stiffness for equilibrium: differences over the frame coordinates
% only, plus the modal stiffness Kc of the elastic units.
    frames = arrayfun(@(B) B.cols(1:6), s.units, 'UniformOutput', false);
    K = s.stiffness(q, u, lam, [], [frames{:}]);
    for B = s.units(logical([s.units.elastic]))
        c = B.cols(7:end);  K(c, c) = K(c, c) + B.Kc;
    end
end
