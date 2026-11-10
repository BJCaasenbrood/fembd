classdef Body < handle
%BODY Linear-elastic (or rigid) solid discretised with Tet4/Tet10 elements.
%   b = fembd.Body(src, mat) creates a body from
%     src  an .stl file (meshed with Gmsh), an ASCII .msh file, or a struct
%          with fields X (n x 3 nodes) and T (ne x 4|10 tets, Gmsh order)
%     mat  material struct with fields D (6x6 Voigt) and rho, see fembd.material
%          (a rigid body only needs rho)
%
%   Options (name-value):
%     Name     body name (default: file name)
%     Scale    factor applied to the coordinates, e.g. 1e-3 for mm -> m
%     Size, Order, Remesh, FeatureAngle   passed to fembd.meshSTL
%     Rigid    true: non-deformable body with 6 dofs at RefPoint; its mesh
%              is only used for mass, inertia, interfaces and plotting
%     RefPoint reference point of a rigid body (default: centroid)
%     Damping  Rayleigh damping [alpha beta] of an elastic body (default
%              [0 0]): D = alpha*M + beta*K. beta = 2*zeta/omega gives the
%              damping ratio zeta at omega, rising linearly with frequency;
%              alpha also damps the body's rigid-body motion (like a drag)
%
%   Rigid bodies without a volume mesh:
%     b = fembd.Body.rigidSurface(file, rho, ...)  STL surface only (no Gmsh)
%     b = fembd.Body.pointMass(m, c, J, ...)       mass properties only
%
%   Degrees of freedom are node-wise interleaved: [ux1 uy1 uz1 ux2 ...].

properties
    Name (1,:) char
    Mat  struct
    X    (:,3) double   % nodes
    T    double         % tetrahedra (ne x 4 or ne x 10)
    F    double         % outward boundary faces (nf x 3 or nf x 6)
    Rigid (1,1) logical = false
    Ref  (1,3) double   % reference point of a rigid body
    MassProps = []      % struct(m, c, J) of a mesh-free rigid body (J about c)
    Damping (1,2) double = [0 0]   % Rayleigh damping [alpha beta]: D = alpha*M + beta*K
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
        end
        if isstruct(src)
            X = src.X;  T = src.T;
        else
            [~, name, ext] = fileparts(src);
            if isempty(opts.Name), opts.Name = name; end
            if strcmpi(ext, '.msh')
                [X, T] = fembd.readMsh(src);
            else
                [X, T] = fembd.meshSTL(src, Size = opts.Size, Order = opts.Order, ...
                    Remesh = opts.Remesh, FeatureAngle = opts.FeatureAngle);
            end
        end
        assert(opts.Rigid || isfield(mat, 'D'), 'fembd:Body', ['%s: an elastic body needs ' ...
            'a material with stiffness D, e.g. fembd.material(E, nu, rho). For a rigid ' ...
            'body set Rigid = true.'], opts.Name);
        assert(~opts.Rigid || ~any(opts.Damping), 'fembd:Body', ['%s: Damping is ' ...
            'material damping of an elastic body; damp rigid bodies with joint ' ...
            'dampers.'], opts.Name);
        b.Name = opts.Name;
        b.Mat  = mat;
        b.X    = opts.Scale*X;
        b.T    = orient(b.X, T);
        b.Rigid = opts.Rigid;
        b.Damping = opts.Damping;
        if isempty(T)
            b.F = zeros(0, 3);  b.Ref = zeros(1, 3);
        else
            b.F = boundary(b.T);  b.Ref = centroid(b.X, b.T);
        end
        if ~isempty(opts.RefPoint), b.Ref = opts.RefPoint; end
    end

    function n = get.nNode(b), n = size(b.X, 1); end
    function n = get.nDof(b),  n = 3*size(b.X, 1); end

    function transform(b, R, t)
    %TRANSFORM Rigidly move the body: x -> R*x + t.
        b.X = b.X*R' + t(:)';
        b.Ref = b.Ref*R' + t(:)';
        if ~isempty(b.MassProps)
            b.MassProps.c = b.MassProps.c*R' + t(:)';
            b.MassProps.J = R*b.MassProps.J*R';
        end
    end

    function setMass(b, m, opts)
    %SETMASS Set the mass properties of a rigid body (e.g. from CAD).
    %   b.setMass(m) sets the mass and scales the inertia (same distribution).
    %   Options:
    %     CenterOfMass  1x3 centre of mass (default: keep the current one)
    %     Inertia       3x3 inertia tensor, or 1x3 principal moments
    %                   (default: current inertia scaled by the mass ratio)
    %     InertiaPoint  point about which Inertia is given (default: the
    %                   centre of mass), e.g. the CAD part origin
    %     Axes          3x3 axes in which Inertia is given (columns, in the
    %                   assembly frame; default eye(3)), e.g. principal axes
    %     LengthScale   factor for lengths in CenterOfMass, InertiaPoint and
    %                   Inertia, e.g. 1e-3 for mm and kg mm^2 (default 1)
    %   Units otherwise SI, in the assembly frame. A meshed rigid body uses
    %   these values instead of its mesh mass.
        arguments
            b
            m (1,1) double {mustBePositive}
            opts.CenterOfMass double = []
            opts.Inertia double = []
            opts.InertiaPoint double = []
            opts.Axes (3,3) double = eye(3)
            opts.LengthScale (1,1) double {mustBePositive} = 1
        end
        assert(b.Rigid, 'fembd:setMass', '%s: setMass is for rigid bodies only.', b.Name);
        p = b.MassProps;
        if isempty(p)                            % meshed rigid body: from its surface
            [p.m, p.c, p.J] = surfaceMass(b.X, b.F(:, 1:3), b.Mat.rho);
        end
        s = opts.LengthScale;  c = p.c;
        if ~isempty(opts.CenterOfMass), c = s*opts.CenterOfMass(:)'; end
        if isempty(opts.Inertia)
            J = p.J*(m/p.m);
        else
            J = opts.Inertia;
            if isvector(J), J = diag(J); end
            J = s^2*opts.Axes*J*opts.Axes';
            if ~isempty(opts.InertiaPoint)       % parallel-axis theorem to the centre of mass
                r = c - s*opts.InertiaPoint(:)';
                J = J - m*(dot(r, r)*eye(3) - r'*r);
            end
        end
        assert(norm(J - J', 'fro') <= 1e-6*norm(J, 'fro'), 'fembd:setMass', ...
            '%s: the inertia tensor must be symmetric.', b.Name);
        J = (J + J')/2;  Jp = sort(eig(J));
        assert(Jp(1) > 0 && Jp(3) <= (Jp(1) + Jp(2))*(1 + 1e-9), 'fembd:setMass', ...
            ['%s: invalid inertia (principal moments %s): they must be positive and ' ...
             'satisfy J3 <= J1 + J2. Check that it is about the centre of mass.'], ...
            b.Name, mat2str(Jp', 4));
        b.MassProps = struct('m', m, 'c', c, 'J', J);
    end

    function M6 = rigidMass(b)
    %RIGIDMASS 6x6 mass matrix of a mesh-free rigid body about its reference
    %   point, for the dofs [ux uy uz rx ry rz].
        p = b.MassProps;  d = p.c - b.Ref;
        S = [0 -d(3) d(2); d(3) 0 -d(1); -d(2) d(1) 0];
        M6 = [p.m*eye(3), -p.m*S; p.m*S, p.J - p.m*S*S];
    end

    function idx = nodes(b, sel)
    %NODES Indices of nodes for which sel(X) is true (X is n x 3), or the
    %   node indices sel themselves if sel is numeric.
        if isnumeric(sel), idx = sel(:); else, idx = find(sel(b.X)); end
    end

    function [K, M] = matrices(b, opts)
    %MATRICES Body stiffness K and mass M (sparse, nDof x nDof); K = 0 for
    %   a rigid body.
    %   M is consistent, or diagonal with Lumped = true (HRZ: element mass
    %   diagonals scaled to preserve the element mass; row-sum lumping
    %   would give negative Tet10 corner masses).
        arguments
            b
            opts.Lumped (1,1) logical = false
            opts.Chunk (1,1) double = 2e4     % elements per vectorised block
        end
        chunk = opts.Chunk;
        [ne, nn] = size(b.T);  nd = 3*nn;  p = 1 + (nn == 10);
        [xk, wk] = tetQuad(2*p - 2);  [~, dNk] = tetShape(xk, nn);
        [xm, wm] = tetQuad(2*p);      [Nm, dNm] = tetShape(xm, nn);
        D = [];  if ~b.Rigid, D = b.Mat.D; end
        K = sparse(b.nDof, b.nDof);  Ms = sparse(b.nNode, b.nNode);
        bad = false(1, ne);                     % elements with det(J) <= 0
        for e0 = 1:chunk:ne
            Te = b.T(e0:min(e0 + chunk - 1, ne), :)';  nc = size(Te, 2);
            Xe = permute(reshape(b.X(Te, :), nn, nc, 3), [1 3 2]);
            Ke = zeros(nd, nd, nc);  Me = zeros(nn, nn, nc);
            for g = 1:numel(wk)*~b.Rigid
                [G, dJ] = gradients(dNk(:,:,g), Xe);
                B = zeros(6, nd, nc);
                B(1, 1:3:end, :) = G(1,:,:);  B(2, 2:3:end, :) = G(2,:,:);
                B(3, 3:3:end, :) = G(3,:,:);
                B(4, 2:3:end, :) = G(3,:,:);  B(4, 3:3:end, :) = G(2,:,:);
                B(5, 1:3:end, :) = G(3,:,:);  B(5, 3:3:end, :) = G(1,:,:);
                B(6, 1:3:end, :) = G(2,:,:);  B(6, 2:3:end, :) = G(1,:,:);
                Ke = Ke + pagemtimes(B, 'transpose', pagemtimes(D, B), 'none').*(wk(g)*dJ);
            end
            for g = 1:numel(wm)
                [~, dJ] = gradients(dNm(:,:,g), Xe);
                bad(e0:e0+nc-1) = bad(e0:e0+nc-1) | dJ(:)' <= 0;
                Me = Me + (Nm(g,:)'*Nm(g,:)).*(wm(g)*dJ);
            end
            dof = reshape(3*permute(Te, [3 1 2]) - [2; 1; 0], nd, nc);
            K  = K  + sparse(repmat(dof, nd, 1), repelem(dof, nd, 1), Ke(:), b.nDof, b.nDof);
            if opts.Lumped
                Me = reshape(Me, nn*nn, nc);  d = Me(1:nn+1:end, :);
                Ms = Ms + sparse(Te, Te, d.*sum(Me, 1)./sum(d, 1), b.nNode, b.nNode);
            else
                Ms = Ms + sparse(repmat(Te, nn, 1), repelem(Te, nn, 1), Me(:), b.nNode, b.nNode);
            end
        end
        if any(bad)
            warning('fembd:Body', '%s: %d element(s) with det(J) <= 0; results are unreliable.', ...
                    b.Name, nnz(bad));
        end
        K = (K + K')/2;
        M = kron((Ms + Ms')/2*b.Mat.rho, speye(3));
        if b.Rigid && ~isempty(b.MassProps)      % mass from setMass, not from the mesh
            M = sparse(b.nDof, b.nDof);
        end
    end
end

methods (Static)
    function b = rigidSurface(file, rho, opts)
    %RIGIDSURFACE Rigid body from a closed STL surface, without volume mesh.
    %   b = fembd.Body.rigidSurface(file, rho) reads the STL triangles, which
    %   are used for ties and plotting; mass, centre of mass and inertia are
    %   computed exactly from the surface (divergence theorem).
    %   Options: Name, Scale, RefPoint (default: centre of mass).
        arguments
            file (1,:) char
            rho (1,1) double
            opts.Name (1,:) char = ''
            opts.Scale (1,1) double = 1
            opts.RefPoint double = []
        end
        tr = stlread(file);
        X = opts.Scale*tr.Points;  F = tr.ConnectivityList;
        [m, c, J, flip] = surfaceMass(X, F, rho);
        if flip, F = F(:, [1 3 2]); end
        if isempty(opts.Name), [~, opts.Name] = fileparts(file); end
        if isempty(opts.RefPoint), opts.RefPoint = c; end
        b = fembd.Body(struct('X', X, 'T', zeros(0, 4)), struct('rho', rho), ...
                     Name = opts.Name, Rigid = true, RefPoint = opts.RefPoint);
        b.F = F;
        b.MassProps = struct('m', m, 'c', c, 'J', J);
    end

    function b = pointMass(m, c, J, opts)
    %POINTMASS Rigid body given by its mass properties only: mass m,
    %   centre of mass c (1x3) and inertia tensor J (3x3, about c; default
    %   0). It has no surface, so connect it with a joint (a.joint).
    %   Options: Name (default 'mass'), RefPoint (default c).
        arguments
            m (1,1) double
            c (1,3) double
            J (3,3) double = zeros(3)
            opts.Name (1,:) char = 'mass'
            opts.RefPoint double = []
        end
        if isempty(opts.RefPoint), opts.RefPoint = c; end
        b = fembd.Body(struct('X', zeros(0, 3), 'T', zeros(0, 4)), struct('rho', 0), ...
                     Name = opts.Name, Rigid = true, RefPoint = opts.RefPoint);
        b.MassProps = struct('m', m, 'c', c, 'J', (J + J')/2);
    end
end
end

function [m, c, J, flip] = surfaceMass(X, F, rho)
% Mass, centre of mass and inertia (about c) of a closed triangulated
% surface, summing signed tetrahedra (origin, a, b, c).
o = mean(X, 1);  P = X - o;
a = P(F(:,1),:);  b = P(F(:,2),:);  e = P(F(:,3),:);  s = a + b + e;
d = dot(a, cross(b, e, 2), 2);                  % 6 x signed volume
V = sum(d)/6;
S1 = sum(d.*s, 1)/24;                           % int x dV
S2 = (a'*(d.*a) + b'*(d.*b) + e'*(d.*e) + s'*(d.*s))/120;   % int x x' dV
flip = V < 0;
if flip, [V, S1, S2] = deal(-V, -S1, -S2); end
m = rho*V;  r = S1/V;  c = o + r;
J = rho*(trace(S2)*eye(3) - S2) - m*(dot(r, r)*eye(3) - r'*r);
end

function [G, dJ] = gradients(dN, Xe)
% Physical shape-function gradients G (3 x nn x ne) and det(J) (1 x 1 x ne).
J  = pagemtimes(dN, 'transpose', Xe, 'none');
dJ = J(1,1,:).*(J(2,2,:).*J(3,3,:) - J(2,3,:).*J(3,2,:)) ...
   - J(1,2,:).*(J(2,1,:).*J(3,3,:) - J(2,3,:).*J(3,1,:)) ...
   + J(1,3,:).*(J(2,1,:).*J(3,2,:) - J(2,2,:).*J(3,1,:));
G  = pagemtimes(pageinv(J), dN');
end

function c = centroid(X, T)
% Volume centroid (straight-sided tets).
a = X(T(:,2),:) - X(T(:,1),:);  b = X(T(:,3),:) - X(T(:,1),:);
v = dot(cross(a, b, 2), X(T(:,4),:) - X(T(:,1),:), 2);
c = v'*(X(T(:,1),:) + X(T(:,2),:) + X(T(:,3),:) + X(T(:,4),:))/4/sum(v);
end

function T = orient(X, T)
% Flip negatively oriented tets.
a = X(T(:,2),:) - X(T(:,1),:);  b = X(T(:,3),:) - X(T(:,1),:);
c = X(T(:,4),:) - X(T(:,1),:);
neg = dot(cross(a, b, 2), c, 2) < 0;
perm = [1 3 2 4 7 6 5 8 10 9];
T(neg,:) = T(neg, perm(1:size(T,2)));
end

function F = boundary(T)
% Outward boundary faces: faces that belong to exactly one tet.
if size(T, 2) == 4
    loc = [1 3 2; 1 2 4; 1 4 3; 2 3 4];
else
    loc = [1 3 2 7 6 5; 1 2 4 5 10 8; 1 4 3 8 9 7; 2 3 4 6 9 10];
end
F = reshape(permute(reshape(T(:, loc'), [], size(loc, 2), 4), [1 3 2]), [], size(loc, 2));
[~, i, j] = unique(sort(F(:, 1:3), 2), 'rows');
F = F(i(accumarray(j, 1) == 1), :);
end
