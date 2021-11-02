function [U,error_all] = ADAM_epoch_unbiased_norm(prec,U,X)

rng(12);

% X = double(tensor(X));

s = size(X);
N = length(s);
r = size(U{1},2);

M = 5*ones(N,1);
% M = ceil(s*0.1);
% M = [10,5,2];

fp.format = 'c';
fp.params = [4,7] ;
fp.round = 5;

% fp.format = 'h';

alpha = 1e-4;
gamma = 1e-5;
beta_1 = 0.9;
beta_2 = 0.999;

eta = 0.8;
epsilon = 1e-8;
wd = 0.001;

error_all = [];
error_w_all = [];
v = cell(N,1);

max_k = 0;

for j = 1:N
   v{j} = zeros(size(U{j})) ;
   m{j} = zeros(s(j),r);
   num{j} = zeros(s(j),1);
end

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

if prec == 0
%     v = cellfun(@(x)half(x),v,'UniformOutput',0);
%     U = cellfun(@(x)half(x),U,'UniformOutput',0);
end

U_w = cellfun(@(x)zeros(size(x)),U,'UniformOutput',0);
t_w = 0;
flag = false;
fail = 0;


for t = 1:300
    
   ind_perm = randperm(D);
   ind_num = 0;
   for j = 1:N
%        v{j} = zeros(size(U{j})) ;
%        m{j} = zeros(s(j),r);
%        num{j} = zeros(s(j),1);
    end
   
   tic,
%    alpha = 0.1/t;
   for batch = ind_perm
       ind_num = ind_num + 1;
       
       c = cell(N,1);
       [c{:}] = ind2sub(n_d,batch);
       
       mul = 1;
       n = cellfun(@(x,y)x{y}, n_all,c,'UniformOutput',false);
       U_tmp = cellfun(@(x,y)y(x,:)/mul, n,U,'UniformOutput',false);
%        U_tmp = cellfun(@(x,y)y(x,:), n,U,'UniformOutput',false);

       X_tmp = X(n{:})/mul^N;
       G = gradient_full_sto_fp(prec,U_tmp,X_tmp,fp);
       G = cellfun(@(x)x*mul^N, G,'UniformOutput',false);
       

       for k = 1:max_k
           G_tmp = gradient_full_sto_fp(prec,U_tmp,X_tmp,fp);
           G = cellfun(@(x,y)(k*x+y)/(k+1), G, G_tmp,'UniformOutput',false);
       end
%        G = gradient_full_sto(prec,U_tmp,X_tmp);
       
       s_tmp = size(X_tmp);
       p = num2cell(prod(s_tmp)./s_tmp.');
       
       for j = 1:N
            num{j}(n{j}) = num{j}(n{j}) + p{j};
           
            m{j}(n{j},:) = beta_1*m{j}(n{j},:) + (1-beta_1)*double(G{j})/p{j};
            v{j}(n{j},:) = beta_2*v{j}(n{j},:) + (1-beta_2)*(double(G{j})/p{j}).^2; 
            
            m_tilde{j}(n{j},:) = m{j}(n{j},:)./(1-beta_1.^num{j}(n{j}));
            v_tilde{j}(n{j},:) = v{j}(n{j},:)./(1-beta_2.^num{j}(n{j}));
            
            tmp = vecnorm(U{j}).^2;
            tmp = tmp(:,1:r-1) - tmp(:,2:r);
            tmp_1 = [tmp.*U{j}(:,1:r-1), zeros(s(j),1)];
            tmp_2 = [zeros(s(j),1),tmp.*U{j}(:,2:r)];
            
%             tmp_1 = [tmp, zeros(s(j),1)];
%             tmp_2 = [zeros(s(j),1), tmp];
            
            U{j}(n{j},:) = U{j}(n{j},:) - alpha*(m_tilde{j}(n{j},:)./(sqrt(v_tilde{j}(n{j},:))+epsilon));
            U{j} = U{j} - gamma * (tmp_1 - tmp_2);
       end

       
      
       
       if flag == true && mod(ind_num,2) == 0
            nU = cellfun(@(x)double(x),U,'UniformOutput',0);
            U_w =  cellfun(@(x,y)(t_w*y+x)/(t_w+1),nU,U_w,'UniformOutput',0);
            t_w = t_w + 1;
       end
       


        
   end
   nU = cellfun(@(x)double(x),U,'UniformOutput',0);
   nX = ktensor(nU);
   
%    U_all{t} = nU;

   error = norm(minus(full(X),full(nX)));
   error_all = [error_all,error];
   
%    alpha = 0.1/t;
   if t>1 && error > error_all(t-1)*1.0
       fail = fail + 1;
       if fail == 1
            fail = 0;
            alpha = max(alpha * 0.1, 0.1^5);
       end
      
   end
%    alpha = 0.1/t;
   
   normX = norm(X(:));
   if error/normX<=1e-3 
%       U_w =  cellfun(@(x,y)(t_w*y+x)/(t_w+1),nU,U_w,'UniformOutput',0);
%       t_w = t_w + 1;
      flag = true;
%       U_w = U;
%       U = U_w;
   end
   
   nX_w = ktensor(U_w);
   error_w = norm(minus(full(X),full(nX_w)));
   error_w_all = [error_w_all, error_w];
   
   if flag == true && error_w > error*1.1
%        U_w = U;
%        t_w = 1;
   end
   
   disp(['epoch = ', num2str(t)]);
   disp(['error = ', num2str(error)]);
   disp(['error_w = ', num2str(error_w)]);
   disp(['time = ', num2str(toc)]);
   disp([' ']);
%    toc,
   
%    t,
%    error, 
   

end


function v = update_v(G,v,p,n)
%     v(n,:) = alpha*G/p + eta*v(n,:);
%     v(n,:) = min(max(alpha*double(G)/p + eta*v(n,:),-1),1);
    v(n,:) = min(max(alpha*(G)/p + eta*v(n,:),-1),1);
%     v(n,:) = min(max(alpha*G + eta*v(n,:),-1),1);
end

function U = update_U(U,v,n)
    U(n,:) = U(n,:) - v(n,:);
end
end

