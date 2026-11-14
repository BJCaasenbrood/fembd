function B = addReshapes(B)
% Reshaped inertia invariants for elasticTerms (computed once).
ne = B.ne;
B.Q9  = reshape(B.Q, 9, []); % Q_ab, (a,b) pairs
B.CRx = reshape(permute(B.CR, [1 3 2]), 3*ne, 1 + ne);
B.Qaz = reshape(B.Q(:,:,:,2:end), 9, []); % Q_{a,j+1}
B.Qju = reshape(B.Q(:,:,2:end,:), 9*ne, 1 + ne); % Q_{j+1,b}
B.Qjz = reshape(B.Q(:,:,2:end,2:end), 9*ne, ne); % Q_{j+1,k+1}
end
