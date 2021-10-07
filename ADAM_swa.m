function [U,error_all,error_swa] = ADAM_swa(prec,U,X)
%ADAM with stochasitc weight averaging
rng(15);

s = size(X);
M = prod(s)/100;

beta_1 = 0.9;
beta_2 = 0.999;
epsilon = 1e-5;
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

U_w = U;

if prec == 0
    U = cellfun(@(x)half(x),U,'UniformOutput',0);
end


w = 0;
error_all = [];
error_swa = [];

for t = 10000:30000
    n = [];
    for j = 1:N
        tmp = randi(s(j),M,1);
        n = [n,tmp];
    end
    G = gradient_M(prec,U,X,n);
    for j = 1:N
       m{j} = beta_1*m{j} + (1-beta_1)*G{j};
       v{j} = beta_2*v{j} + (1-beta_2)*G{j}.^2;
       
       m_tilde{j} = m{j}/(1-beta_1^t);
       v_tilde{j} = v{j}/(1-beta_2^t);
       
       U{j} = U{j} - alpha*m_tilde{j}./(sqrt(v_tilde{j})+epsilon);
    end
    nU = cellfun(@(x)double(x),U,'UniformOutput',0);
    nX = ktensor(nU);

    error = norm(minus(full(X),full(nX)));
    error_all = [error_all,error];
        
    if t>10100 && error>1.01*error_all(t-10001)
       alpha = 0.2 * alpha; 
    end
    
%     if mod(t,100)==0 && t>8000
%         for j = 1:N
%            U_w{j} = (U_w{j}*w+nU{j})/(w+1);
%            w = w + 1;
%            nX = ktensor(U_w);
%            error = norm(minus(full(X),full(nX)));
%            error_swa = [error_swa,error];
%         end
%     end
end