function Xn = n_mode_unfold(X,n)
N = ndims(X);
sizes = size(X);
Xn = permute(X,[n,1:n-1,n+1:N]);
Xn = reshape(Xn, sizes(n), numel(X)/sizes(n));