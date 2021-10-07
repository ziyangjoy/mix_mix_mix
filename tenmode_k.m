function M = tenmode_k(X,k)
s = size(X);
N = length(s);

X = permute(X, [k,1:k-1,k+1:N]);
M = X(:,:).';


