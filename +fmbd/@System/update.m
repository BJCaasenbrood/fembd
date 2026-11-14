function update(s)
%UPDATE Rebuild the system on its next use, after changing a body
%   directly (e.g. b.Modes, b.Damping, b.setMass); changes made through
%   the System methods are picked up automatically.
    s.built = false;
end
