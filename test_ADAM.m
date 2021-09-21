%load aminoacids.mat
clear all;

rng(12);

N = 3;

prec = 0;

size = 20;
r = 20;

s = [size/2,size*2,size];

A = cell(N,1);
for i = 1:N
    A{i} = (randi(11,s(i),r)-6);
%     A{i} = double(A{i});
end
X = ktensor(A);
X = tensor(double(X)/max(abs(double(X(:)))));

iterCG = 10;
iterSG = 50;
U = cell(N,1);


for i = 1:N
    U{i} = randn(s(i),r);
end

% if prec == 0
%     U = cellfun(@(x)half(x),U,'UniformOutput',0);
% end


D = prod(s);
M = 30;

beta_1 = 0.9;
beta_2 = 0.999;
epsilon = 1e-5;
alpha = 0.01;

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
       m{j} = beta_1*m{j} + (1-beta_1)*G{j}/M;
       v{j} = beta_2*v{j} + (sqrt(1-beta_2)*G{j}/M).^2;
       
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

U = cellfun(@(x)double(x),U,'UniformOutput',0);
nX = ktensor(U);

error = norm(minus(full(nX),full(X)));