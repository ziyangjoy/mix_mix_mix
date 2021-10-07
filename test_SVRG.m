%load aminoacids.mat
clear all;

addpath('./tensor_toolbox-v3.2');

rng(15);

N = 3;

size = 70;
r = 20;

s = [size/2,size*2,size];

% s = [200 150 100 50];
% r = 5;

A = cell(N,1);
for i = 1:N
%     A{i} = double(half(randn(s(i),r)));
    A{i} = randn(s(i),r)+.1;
end
X = ktensor(A);
X = full(X);
% X = X/max(X(:));



U = cell(N,1);

rng(22);

for i = 1:N
    U{i} = randn(s(i),r);
end

D = prod(s);
prec = 2;
M = floor(D/500);
M = 200;

alpha = 0.005;



error_all = [];

% load('SVRG_test_data.mat')
% alpha = 0.0001;
for i = 1:1000
    n = [];
    for j = 1:N
        tmp = randi(s(j),50*M,1);
        n = [n,tmp];
    end
%     G_full = gradient_full(prec,U,X);
    G_full = gradient_M_new_w(prec,U,X,n);
%     G_full = gradient_full_sto(prec,U,X);
    for j = 1:N
%        G_full{j} = G_full{j}/D;
%        G_full{j} = G_full{j}/D*s(j);
    end
    nU = U;

    for k = 1:50
        n = [];
        for j = 1:N
            tmp = randi(s(j),M,1);
            n = [n,tmp];
        end
%         G = gradient_M(prec,nU,X,n);
%         G_tilde = gradient_M(prec,U,X,n);
        G = gradient_M_new_w(prec,nU,X,n);
        G_tilde = gradient_M_new_w(prec,U,X,n);
%         G = gradient_sto_new_w(prec,nU,X,n);
%         G_tilde = gradient_sto_new_w(prec,U,X,n);
        
        for j = 1:N
%             nU{j} = nU{j} - alpha*(G{j}-G_tilde{j}+G_full{j});
            nU{j} = nU{j} - alpha*(double(G{j})-double(G_tilde{j})+double(G_full{j}));
        end
    end
    U = nU;
    U = cellfun(@(x)double(x),U,'UniformOutput',0);
    nX = ktensor(U);

    error = norm(minus(full(nX),full(X)));
    error_all = [error_all,error];
    
    if mod(i,10)==0
       error, 
    end
end

U = cellfun(@(x)double(x),U,'UniformOutput',0);
nX = ktensor(U);

error = norm(minus(full(nX),full(X)));