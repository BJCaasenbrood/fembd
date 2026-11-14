function [G, dJ] = gradients(dN, Xe)
% Physical shape-function gradients G (3 x nn x ne) and det(J) (1 x 1 x ne).
J  = pagemtimes(dN, 'transpose', Xe, 'none');
dJ = J(1,1,:).*(J(2,2,:).*J(3,3,:) - J(2,3,:).*J(3,2,:)) ...
   - J(1,2,:).*(J(2,1,:).*J(3,3,:) - J(2,3,:).*J(3,1,:)) ...
   + J(1,3,:).*(J(2,1,:).*J(3,2,:) - J(2,2,:).*J(3,1,:));
G  = pagemtimes(pageinv(J), dN');
end
