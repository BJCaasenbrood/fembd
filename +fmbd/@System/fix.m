function fix(s, id, varargin)
%FIX Support:
%   s.fix(ref, dirs)         holds a reference point (rigid body or
%                            coupling) in dirs (default 1:6, [ux uy uz
%                            rx ry rz], world axes)
%   s.fix(joint, dirs)       holds joint coordinates at their Position
%   s.fix(body, sel, dirs)   clamps nodes of an elastic body (dirs
%                            default 1:3); its elastic group is grounded,
%                            so the supports must hold it (e.g. a face)
    if numel(varargin) >= 2 || (isscalar(varargin) && isa(varargin{1}, 'function_handle'))
        sel = varargin{1};  dirs = 1:3;  if numel(varargin) >= 2, dirs = varargin{2}; end
        [kind, k] = s.lookup(id, false);
        assert(strcmp(kind, 'body') && ~s.Bodies{k}.Rigid, 'fmbd:fix', ...
            'Node supports are for elastic bodies; support a rigid body with s.fix(body, dirs).');
        nodes = s.Bodies{k}.nodes(sel);
        assert(~isempty(nodes), 'fmbd:fix', 'No nodes selected.');
        s.Fixes(end+1) = struct('kind', 'nodes', 'idx', k, 'dirs', dirs, 'nodes', nodes);
    else
        dirs = [];  if isscalar(varargin), dirs = varargin{1}; end
        [kind, k] = s.lookup(id, true);
        assert(ismember(kind, {'joint', 'ref'}), 'fmbd:fix', ['%s is not a joint or reference ' ...
            'point; clamp nodes of an elastic body with s.fix(body, sel).'], string(id));
        s.Fixes(end+1) = struct('kind', kind, 'idx', k, 'dirs', dirs, 'nodes', []);
    end
    s.built = false;
end
