function output(s, name, id, sel, dirs, opts)
%OUTPUT Add a named output.
%   s.output(name, id, sel, dirs) measures joint coordinates, reference-point
%   coordinates or the mean displacement of a node set, selected as in input.
%   Option Type: 'displacement' (default), 'velocity' or 'acceleration'.
    arguments
        s
        name (1,:) char
        id
        sel = []
        dirs = []
        opts.Type (1,:) char {mustBeMember(opts.Type, {'displacement', 'velocity', 'acceleration'})} = 'displacement'
    end
    s.Outputs = [s.Outputs, signalSpec(s, name, id, sel, dirs, opts.Type)];
    s.built = false;
end
