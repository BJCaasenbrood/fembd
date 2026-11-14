classdef Body < handle
%BODY Rigid or linear-elastic body.
%   b = fmbd.Body.elastic(file, mat, ...)   elastic, meshed with Gmsh
%   b = fmbd.Body.rigid(file, rho, ...)     rigid, from its STL surface
%   b = fmbd.Body.pointMass(m, c, J, ...)   rigid, from mass properties
%
%   b = fmbd.Body(src, mat) takes an .stl, an ASCII .msh, or a struct with X
%   (nodes) and T (Tet4/Tet10 in Gmsh order). mat holds D and rho.
%   Options: Name, Scale (e.g. 1e-3 for mm), Size, Order, Remesh,
%   FeatureAngle (meshing), Rigid, RefPoint, Modes (default 10) and
%   Damping = [alpha beta] (D = alpha*M + beta*K).

properties
    Name (1,:) char
    Mat  struct
    X    (:,3) double % nodes
    T    double % tetrahedra (ne x 4 or ne x 10)
    F    double % outward boundary faces (nf x 3 or nf x 6)
    Rigid (1,1) logical = false
    Ref  (1,3) double % reference point of a rigid body
    MassProps = [] % struct(m, c, J) of a mesh-free rigid body (J about c)
    Damping (1,2) double = [0 0] % Rayleigh damping [alpha beta]: D = alpha*M + beta*K
    Modes (1,1) double = 10 % fixed-interface modes of an elastic body
end

properties (Dependent)
    nNode
    nDof
end

methods
    function b = Body(src, mat, opts)
        arguments
            src
            mat struct
            opts.Name (1,:) char = ''
            opts.Scale (1,1) double = 1
            opts.Size (1,1) double = inf
            opts.Order (1,1) double = 2
            opts.Remesh (1,1) logical = true
            opts.FeatureAngle (1,1) double = 40
            opts.Rigid (1,1) logical = false
            opts.RefPoint double = []
            opts.Damping (1,2) double {mustBeNonnegative} = [0 0]
            opts.Modes (1,1) double {mustBeInteger, mustBeNonnegative} = 10
        end
        if isstruct(src)
            X = src.X;  T = src.T;
        else
            [~, name, ext] = fileparts(src);
            if isempty(opts.Name), opts.Name = name; end
            if strcmpi(ext, '.msh')
                [X, T] = fmbd.readMsh(src);
            else
                [X, T] = fmbd.meshSTL(src, Size = opts.Size, Order = opts.Order, ...
                    Remesh = opts.Remesh, FeatureAngle = opts.FeatureAngle);
            end
        end
        assert(opts.Rigid || isfield(mat, 'D'), 'fmbd:Body', ['%s: an elastic body needs ' ...
            'a material with stiffness D, e.g. fmbd.material(E, nu, rho). For a rigid ' ...
            'body set Rigid = true.'], opts.Name);
        assert(~opts.Rigid || ~any(opts.Damping), 'fmbd:Body', ['%s: Damping is ' ...
            'material damping of an elastic body; damp rigid bodies with joint ' ...
            'dampers.'], opts.Name);
        b.Name = opts.Name;
        b.Mat  = mat;
        b.X    = opts.Scale*X;
        b.T    = orient(b.X, T);
        b.Rigid = opts.Rigid;
        b.Damping = opts.Damping;  b.Modes = opts.Modes;
        if isempty(T)
            b.F = zeros(0, 3);  b.Ref = zeros(1, 3);
        else
            b.F = boundary(b.T);  b.Ref = centroid(b.X, b.T);
        end
        if ~isempty(opts.RefPoint), b.Ref = opts.RefPoint; end
    end

    function n = get.nNode(b), n = size(b.X, 1); end
    function n = get.nDof(b),  n = 3*size(b.X, 1); end
end

methods
    transform(b, R, t)
    setMass(b, m, opts)
    [m, c, J] = massProperties(b)
    idx = nodes(b, sel)
    [K, M] = matrices(b, opts)
end

methods (Static)
    b = elastic(file, mat, varargin)
    b = rigid(file, rho, opts)
    b = pointMass(m, c, J, opts)
end
end
