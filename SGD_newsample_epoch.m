function [U,error_all] = SGD_newsample_epoch(prec,U,X)

rng(12);

s = size(X);
N = length(s);
r = size(U{1},2);

% M = 5*ones(N,1);
M = ceil(s*0.1);


alpha = 0.01;




error_all = [];
v = cell(N,1);
eta = 0.8;



n_all = cell(N,1);
n_d = [];
for j = 1:N
    tmp = randperm(s(j));
    i_st = 1;
    k = 0;
    while i_st<=s(j)
        k = k + 1;
        i_ed = min(i_st + M(j)-1,s(j));
        n_all{j}{k} = tmp(i_st:i_ed);
        i_st = i_ed + 1;
    end
    n_d = [n_d, k];
end

D = prod(n_d);

for j = 1:N
   v{j} = zeros(size(U{j})) ;
end

if prec == 0
    v = cellfun(@(x)half(x),v,'UniformOutput',0);
    U = cellfun(@(x)half(x),U,'UniformOutput',0);
end



for t = 1:300
    
    ind_perm = randperm(D);

   for batch = ind_perm
       c = cell(N,1);
       [c{:}] = ind2sub(n_d,batch);
       n = cell(N,1);
       for j = 1:N
            n{j} = n_all{j}{c{j}};
            U_tmp{j} = U{j}(n{j},:);
       end
       X_tmp = X(n{:});
       G = gradient_full(prec,U_tmp,X_tmp);
       s_tmp = size(X_tmp);
        for j = 1:N
           p = prod(s_tmp([1:j-1,j+1:N]));
%            p = 1;
           G_tmp = G{j};
%            G_tmp = min(max(G{j},-10),10);
           v{j}(n{j},:) = alpha*(G_tmp)/p + eta*v{j}(n{j},:);
           v{j} = min(max(v{j},-1),1);
           U{j}(n{j},:) = U{j}(n{j},:) - v{j}(n{j},:);
        end
        
        nU = cellfun(@(x)double(x),U,'UniformOutput',0);
        nX = ktensor(nU);

        error = norm(minus(full(X),full(nX)));
        error_all = [error_all,error];
        
   end

    
   t,
   error, 
   
%    if error/norm(X)<1e-5
%       break; 
%    end
end



