function input(s, name, id, sel, dirs)
%INPUT Add a named input.
%   s.input(name, joint)            force or torque on the joint coordinates
%   s.input(name, ref, [], dirs)    force (1:3) or moment (4:6) on a reference point
%   s.input(name, body, sel, dirs)  total force on a node set
    if nargin < 4, sel = []; end
    if nargin < 5, dirs = []; end
    s.Inputs = [s.Inputs, s.signalSpec(name, id, sel, dirs, 'force')];
    s.built = false;
end
