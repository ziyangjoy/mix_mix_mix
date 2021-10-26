function [U,error_all] = SVRG_epoch(prec,U,X)

rng(12);

% X = double(tensor(X));

s = size(X);
N = length(s);
r = size(U{1},2);

% M = 5*ones(N,1);
M = ceil(s*0.1);
% M = [10,5,2];



alpha = 0.02;




error_all = [];
v = cell(N,1);

eta = 0.0;



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

U_w = cellfun(@(x)zeros(size(x)),U,'UniformOutput',0);
t_w = 0;
flag = false;
U_old = U;

for t = 1:300
    
   ind_perm = randperm(D);
   ind_num = 0;
   
   tic,
   for batch = ind_perm
       ind_num = ind_num + 1;
       
       if mod(ind_num,100) == 1
           G_full = gradient_full(prec,U,X);
           U_old = U;
       end
       
       c = cell(N,1);
       [c{:}] = ind2sub(n_d,batch);
       n = cellfun(@(x,y)x{y}, n_all,c,'UniformOutput',false);
       U_tmp = cellfun(@(x,y)y(x,:), n,U,'UniformOutput',false);
       U_tmp_old = cellfun(@(x,y)y(x,:), n,U_old,'UniformOutput',false);
%        U_tmp = cellfun(@(x,y)y(x,:), n,U,'UniformOutput',false);

       X_tmp = X(n{:});
       G = gradient_full(prec,U_tmp,X_tmp);
       G_tilde = gradient_full(prec,U_tmp_old,X_tmp);
%        G = gradient_full(prec,U_tmp,X_tmp);
       s_tmp = size(X_tmp);
       p = num2cell(prod(s_tmp)./s_tmp.');
       
       v = cellfun(@update_v, G,v,p,n,G_tilde,G_full,'UniformOutput',false);
       U = cellfun(@update_U,U,v,n,'UniformOutput',false);
       
       if flag == true && mod(ind_num,2) == 0
            nU = cellfun(@(x)double(x),U,'UniformOutput',0);
            U_w =  cellfun(@(x,y)(t_w*y+x)/(t_w+1),nU,U_w,'UniformOutput',0);
            t_w = t_w + 1;
       end
       
%        if mod(ind_num,100) == 0
%            
%            nU = cellfun(@(x)double(x),U,'UniformOutput',0);
%            nX = ktensor(nU);
%            error = norm(minus(full(X),full(nX)));
%            
% %            nX = gen_ten(U);
% %            error = norm(double(X(:)-nX(:)));
%            error_all = [error_all,error];
%        end


        
   end
   nU = cellfun(@(x)double(x),U,'UniformOutput',0);
   nX = ktensor(nU);
   
%    U_all{t} = nU;

   error = norm(minus(full(X),full(nX)));
   error_all = [error_all,error];
   
   normX = norm(X(:));
   if error/normX<=1e-3
%       U_w =  cellfun(@(x,y)(t_w*y+x)/(t_w+1),nU,U_w,'UniformOutput',0);
%       t_w = t_w + 1;
      flag = true;
%       U = U_w;
   end
   
   nX_w = ktensor(U_w);
   error_w = norm(minus(full(X),full(nX_w)));
   
   disp(['epoch = ', num2str(t)]);
   disp(['error = ', num2str(error)]);
   disp(['error_w = ', num2str(error_w)]);
   disp(['time = ', num2str(toc)]);
   disp([' ']);
%    toc,
   
%    t,
%    error, 
   

end


function v = update_v(G,v,p,n,G_tilde,G_full)
    v(n,:) = alpha*(G/p + G_tilde/p - G_full(n,:)/2500) + eta*v(n,:);
%     v(n,:) = min(max(alpha*double(G)/p + eta*v(n,:),-1),1);
    v(n,:) = min(max(v(n,:),-1),1);
%     v(n,:) = min(max(alpha*G + eta*v(n,:),-1),1);
end

function U = update_U(U,v,n)
    U(n,:) = U(n,:) - v(n,:);
end
end

