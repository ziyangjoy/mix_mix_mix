function [U,error_all] = ADAM(prec,U,X)

rng(15);

s = size(X);
M = floor(prod(s)/500);
M = 100;

beta_1 = 0.9;
beta_2 = 0.999;
epsilon = 1e-8;
alpha = .01;

N = length(s);
r = size(U{1},2);
m = cell(N,1);
v = cell(N,1);
m_tilde = cell(N,1);
v_tilde = cell(N,1);

for j = 1:N
    m{j} = zeros(s(j),r);
    v{j} = zeros(s(j),r);
end

error_all = [];


for t = 1:20000
    n = [];
    for j = 1:N
        tmp = randi(s(j),M,1);
        n = [n,tmp];
    end
    G = gradient_M_new_w(prec,U,X,n);
%     G = gradient_M(prec,U,X,n);
    for j = 1:N
%        m{j} = beta_1*m{j} + (1-beta_1)*(G{j});
%        v{j} = beta_2*v{j} + (sqrt(1-beta_2)*(G{j})).^2;
%        m{j} = beta_1*m{j} + (1-beta_1)*double(G{j});
%        v{j} = beta_2*v{j} + (1-beta_2)*(double(G{j})).^2;
       
%        m_tilde{j} = m{j}/(1-beta_1^t);
%        v_tilde{j} = v{j}/(1-beta_2^t);
%        
%        U{j} = U{j} - alpha*m_tilde{j}./(sqrt(v_tilde{j})+epsilon);
       U{j} = U{j} - alpha*(G{j});
    end
    nU = cellfun(@(x)double(x),U,'UniformOutput',0);
    nX = ktensor(nU);

    error = norm(minus(full(X),full(nX)));
    error_all = [error_all,error];
    if mod(t,500)==0
       error, 
    end
end