function v = cr(Q)
% sum_pq Ms_pq y_p x z_q from Q = sum_pq Ms_pq y_p z_q'.
v = [Q(2,3) - Q(3,2); Q(3,1) - Q(1,3); Q(1,2) - Q(2,1)];
end
