function n = names(s)
% All names in use (bodies, reference points, joints).
    n = [cellfun(@(b) b.Name, s.Bodies, 'UniformOutput', false), {s.Refs.Name}, {s.Joints.Name}];
end
