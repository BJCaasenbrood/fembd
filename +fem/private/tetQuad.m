function [xi, w] = tetQuad(p)
%TETQUAD Quadrature on the reference tet, exact for polynomials of degree p.
%   Collapsed (Duffy) Gauss-Legendre rule; weights sum to 1/6.

if p <= 1
    w = 1/6;
    xi = [1 1 1]/4;  
    return
end
n = ceil((p + 3)/2);

% Golub-Welsch
b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);              
[V, L] = eig(diag(b, 1) + diag(b, -1));
x = (diag(L) + 1)/2;  wx = V(1,:)'.^2;
[u, v, s] = ndgrid(x);
w  = reshape(wx.*wx'.*reshape(wx, 1, 1, []), [], 1);
xi = [u(:), v(:).*(1 - u(:)), s(:).*(1 - u(:)).*(1 - v(:))];
w  = w.*(1 - u(:)).^2.*(1 - v(:));
end
