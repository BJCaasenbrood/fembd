function K = frameStiffness(s, q, u, lam)
% Approximate tangent stiffness for equilibrium: differences only over
% the frame coordinates (6 per unit), and the modal stiffness Kc of the
% elastic units exactly. It leaves out how gravity and constraint
% forces change with the modes, which only slows Newton slightly.
    frames = arrayfun(@(B) B.cols(1:6), s.units, 'UniformOutput', false);
    K = s.stiffness(q, u, lam, [], [frames{:}]);
    for B = s.units(logical([s.units.elastic]))
        c = B.cols(7:end);  K(c, c) = K(c, c) + B.Kc;
    end
end
