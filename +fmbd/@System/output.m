function output(s, name, id, sel, dirs, opts)
%OUTPUT Named output (sensor), one signal per direction: joint
%   coordinates, reference-point coordinates or the average displacement
%   of a node set, selected as in input. Type: 'displacement' (default),
%   'velocity' or 'acceleration'.
    arguments
        s
        name (1,:) char
        id
        sel = []
        dirs = []
        opts.Type (1,:) char {mustBeMember(opts.Type, {'displacement', 'velocity', 'acceleration'})} = 'displacement'
    end
    s.Outputs = [s.Outputs, s.signalSpec(name, id, sel, dirs, opts.Type)];
    s.built = false;
end
