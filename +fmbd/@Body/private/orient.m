function T = orient(X, T)
% Flip negatively oriented tets.
a = X(T(:,2),:) - X(T(:,1),:);  b = X(T(:,3),:) - X(T(:,1),:);
c = X(T(:,4),:) - X(T(:,1),:);
neg = dot(cross(a, b, 2), c, 2) < 0;
perm = [1 3 2 4 7 6 5 8 10 9];
T(neg,:) = T(neg, perm(1:size(T,2)));
end
