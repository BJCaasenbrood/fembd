function pairs = rotationPairs(free)
% Axis pairs (i on A, j on B) of the rotation constraints e_i^A . e_j^B = 0:
% none free: all three; one free axis k: k on B against the other two on
% A (axis aligned); two free: the Cardan pair before the locked axis.
fr = free(4:6);
switch nnz(fr)
    case 0, pairs = [2 3; 3 1; 1 2];
    case 1, k = find(fr);  i = setdiff(1:3, k);  pairs = [i(:), [k; k]];
    case 2, k = find(~fr);  pairs = circshift([1 2 3], -k);  pairs = pairs(1:2);
    otherwise, pairs = zeros(0, 2);
end
end
