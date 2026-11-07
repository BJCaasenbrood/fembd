%% Dofs, Supports, Loads and the System
%
%% Dof layout
% The full displacement vector u contains, in order:
%
% # the node dofs [ux1 uy1 uz1 ux2 ...] of every body, in the order the
%   bodies were added;
% # 6 dofs [ux uy uz rx ry rz] per reference point (rigid bodies and
%   couplings), in creation order;
% # the released coordinates of every joint, in creation order.
%
% |a.dof| returns global dof indices:
%
%   d = a.dof(body, sel, dirs)    % nodes selected by sel, dirs in 1:3
%   d = a.dof(ref, [], dirs)      % reference point, dirs in 1:6
%   d = a.dof(joint)              % joint coordinates
%
% A selector |sel| is a function of the node coordinates, e.g.
% |@(X) X(:,1) < 1e-9|, or a list of node indices.
%
%% Supports
%
%   a.fix(body, sel, dirs)        % nodes (dirs default 1:3)
%   a.fix(ref, [], dirs)          % rigid body or coupling, dirs in 1:6
%   a.fix(joint, [])              % lock a joint
%
%% The system
%
%   [K, M, T] = a.system(Lumped = false)
%
% All constraints are eliminated with u = T*q: tie slaves, nodes of rigid
% bodies and couplings, B sides of joints, cut joints and supports. The
% returned K = T'*Kfull*T and M = T'*Mfull*T are symmetric and sparse.
% |a.assemble()| gives the unconstrained Kfull and Mfull, and
% |a.constraints()| gives T and the independent dofs.
%
% |a.qdof(d)| converts global dofs into positions in q, so
% u(d) = q(a.qdof(d)).
%
%% Loads and outputs
% Build loads and outputs in the full space and project them:
%
%   f  = zeros(a.nDof, 1);  f(d) = F;      % nodal forces
%   fq = T'*f;                             % reduced load
%   q  = K\fq;  u = T*q;                   % static solution
%   C  = zeros(1, a.nDof);  C(d(1)) = 1;  Cq = C*T;   % output
%
% Gravity is |T'*(Mfull*g)| with g = -9.81 in every z node dof. Support
% reactions are |Kfull*u - f| at the supported dofs. See
% <tut05_loadsStatics.html Tutorial 5>.
%
%% Plotting
% |a.plot(u, scale)| draws all bodies deformed by the full vector u.
% Since u also holds rotation dofs (radians), scale with the node dofs
% only: |scale = s/max(abs(u(1:3*size(a.X, 1))))|.
