function pairs = rotationPairs(free)
% Axis pairs of the rotation constraints e_i(A) . e_j(B) = 0 of a joint.
fr = free(4:6);
switch nnz(fr)
    case 0, pairs = [2 3; 3 1; 1 2];
    case 1, k = find(fr);  i = setdiff(1:3, k);  pairs = [i(:), [k; k]];
    case 2, k = find(~fr);  pairs = circshift([1 2 3], -k);  pairs = pairs(1:2);
    otherwise, pairs = zeros(0, 2);
end
end
