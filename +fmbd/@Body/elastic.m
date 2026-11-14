function b = elastic(file, mat, varargin)
%ELASTIC Elastic body meshed with Gmsh.
%   b = fmbd.Body.elastic(file, mat, Size = h, Modes = 10, ...). See fmbd.Body.
    b = fmbd.Body(file, mat, varargin{:});
    assert(~b.Rigid, 'fmbd:Body', 'Use fmbd.Body.rigid for rigid bodies.');
end
