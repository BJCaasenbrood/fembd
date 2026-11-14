function update(s)
%UPDATE Rebuild the system after changing a body directly.
%   s.update() is needed after e.g. b.Modes = 20 or b.setMass(...).
    s.built = false;
end
