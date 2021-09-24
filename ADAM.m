function [U,error_all] = ADAM(prec,U,X)

rng(15);

s = size(X);
M = 50;

beta_1 = 0.9;
beta_2 = 0.999;
epsilon = 1e-5;
alpha = 0.01;

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
    G = gradient_M(prec,U,X,n);
    for j = 1:N
       m{j} = beta_1*m{j} + (1-beta_1)*double(G{j});
       v{j} = beta_2*v{j} + (1-beta_2)*(double(G{j})).^2;
       
       m_tilde{j} = m{j}/(1-beta_1^t);
       v_tilde{j} = v{j}/(1-beta_2^t);
       
       U{j} = U{j} - alpha*m_tilde{j}./(sqrt(v_tilde{j})+epsilon);
    end
    nU = cellfun(@(x)double(x),U,'UniformOutput',0);
    nX = ktensor(nU);

    error = norm(minus(full(X),full(nX)));
    if error<1e-2
        break;
    end
    error_all = [error_all,error];
end