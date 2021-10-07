function [U,error_all] = SGD_epoch_GPU(prec,U,X)

% rng(12);

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

v = cellfun(@(x)gpuArray(x),v,'UniformOutput',0);
U = cellfun(@(x)gpuArray(x),U,'UniformOutput',0);

X = gpuArray(double(tensor(X)));

for t = 1:1
    
   ind_perm = randperm(D);
   
   ind_num = 0;
   tic,
   for batch = ind_perm
       ind_num = ind_num + 1;
       
       c = cell(N,1);
       [c{:}] = ind2sub(n_d,batch);
       n = cellfun(@(x,y)x{y}, n_all,c,'UniformOutput',false);
       U_tmp = cellfun(@(x,y)y(x,:), n,U,'UniformOutput',false);

       X_tmp = X(n{:});
       G = gradient_full_ownfunc(prec,U_tmp,X_tmp);
       s_tmp = size(X_tmp);
       p = num2cell(prod(s_tmp)./s_tmp.');


       v = cellfun(@update_v, G,v,p,n,'UniformOutput',false);
       U = cellfun(@update_U,U,v,n,'UniformOutput',false);
       
       if mod(ind_num,1000) == 0
           nU = cellfun(@(x)gather(double(x)),U,'UniformOutput',0);
           nX = gpuArray(double(tensor(ktensor(nU))));

%            error = norm(minus(full(X),full(nX)));
           error = norm(X(:)-nX(:));
           error_all = [error_all,error];
       end

%        nU = cellfun(@(x)double(x),U,'UniformOutput',0);
%        nX = ktensor(nU);
% 
%        error = norm(minus(full(X),full(nX)));
%        error_all = [error_all,error];
        
   end
   toc,
   t,
   error, 
   
%    if error/norm(X)<1e-5
%       break; 
%    end
end


function v = update_v(G,v,p,n)
%     v = eta*v;
    v(n,:) = alpha*G/p + eta*v(n,:);
end

function U = update_U(U,v,n)
    U(n,:) = U(n,:) - v(n,:);
end
end

