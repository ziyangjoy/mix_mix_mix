%load aminoacids.mat
clear all;

rng(15);

N = 3;

size = 20;
r = 10;

s = [size/2,size*2,size];

A = cell(N,1);
for i = 1:N
%     A{i} = randn(s(i),r);
    A{i} = (randi(11,s(i),r)-6);
end
X = ktensor(A);
X = tensor(double(X)/max(abs(double(X(:)))));
% load('fail_X.mat')

iterCG = 10;
iterSG = 50;
U = cell(N,1);

rng(22);
for i = 1:N
    U{i} = randn(s(i),r);
end

D = prod(s);
prec = 0;
M = 10;
alpha = 0.1;

lambda = 0.001;

error_all = [];
for i = 1:500
    G_full = gradient_full(prec,U,X);
    for j = 1:N
       G_full{j} = G_full{j}/D;
    end
    nU = U;
%     alpha = alpha * 0.999;
    for k = 1:20
        n = [];
        for j = 1:N
            tmp = randi(s(j),M,1);
            n = [n,tmp];
        end
        G = gradient_M(prec,nU,X,n);
        G_tilde = gradient_M(prec,U,X,n);
        
        for j = 1:N
%             nU{j} = nU{j} - alpha*(G{j}/M);
            nU{j} = nU{j} - alpha*(G{j}/M-G_tilde{j}/M+G_full{j});
        end
    end
    U = nU;
    U = cellfun(@(x)double(x),U,'UniformOutput',0);
    nX = ktensor(U);

    error = norm(minus(full(nX),full(X)));
    error_all = [error_all,error];
end

U = cellfun(@(x)double(x),U,'UniformOutput',0);
nX = ktensor(U);

error = norm(minus(full(nX),full(X)));