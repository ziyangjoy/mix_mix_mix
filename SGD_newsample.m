function [U,error_all] = SGD_newsample(prec,U,X)

% rng(12);

s = size(X);
N = length(s);
r = size(U{1},2);

% M = 5*ones(N,1);
M = ceil(s*0.2);


alpha = .01;




error_all = [];
v = cell(N,1);
eta = 0.9;

for j = 1:N
   v{j} = zeros(size(U{j})) ;
end

for t = 1:50000
    n = cell(N,1);
    for j = 1:N
        n{j} = randsample(s(j),M(j));
        U_tmp{j} = U{j}(n{j},:)-v{j}(n{j},:);
%         U_tmp{j} = U{j}(n{j},:);
    end
    X_tmp = X(n{:});
    G = gradient_full(prec,U_tmp,X_tmp);


    for j = 1:N
%        G{j} = min(max(G{j},-10),10);
       c = prod(M([1:j-1,j+1:N]));
%        c = prod(M);
       v{j}(n{j},:) = alpha*(G{j})/c + eta*v{j}(n{j},:);
       U{j}(n{j},:) = U{j}(n{j},:) - v{j}(n{j},:);
    end
    nU = cellfun(@(x)double(x),U,'UniformOutput',0);
    nX = ktensor(nU);

    error = norm(minus(full(X),full(nX)));
    error_all = [error_all,error];
    
    if mod(t,500)==0
       error, 
    end
    

    
end