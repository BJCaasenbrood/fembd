function input(s, name, id, sel, dirs)
%INPUT Named input, one signal per direction:
%   s.input(name, joint)             force/torque on the joint coordinates
%   s.input(name, ref, [], dirs)     force (1:3) or moment (4:6, world
%                                    axes) on a rigid body or coupling
%   s.input(name, body, sel, dirs)   total force on a node set, split
%                                    equally over the nodes
    if nargin < 4, sel = []; end
    if nargin < 5, dirs = []; end
    s.Inputs = [s.Inputs, s.signalSpec(name, id, sel, dirs, 'force')];
    s.built = false;
end
