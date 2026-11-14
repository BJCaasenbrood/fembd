function b = elastic(file, mat, varargin)
%ELASTIC Elastic body meshed with Gmsh: b = fmbd.Body.elastic(file, mat,
%   Size = ..., Order = 2, Modes = 10, Damping = [alpha beta], ...). See
%   fmbd.Body for the options.
    b = fmbd.Body(file, mat, varargin{:});
    assert(~b.Rigid, 'fmbd:Body', 'Use fmbd.Body.rigid for rigid bodies.');
end
