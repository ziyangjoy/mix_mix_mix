%load aminoacids.mat
clear all;

addpath('./chop');
addpath('./tensor_toolbox-v3.2');

rng(12);

N = 3;

prec = 2;

size = 20;
r = 10;

s = [size/2,size*2,size];

A = cell(N,1);
for i = 1:N
    A{i} = double(half(randn(s(i),r)));
end
X = ktensor(A);
% X = tensor(double(X)/max(abs(double(X(:)))));


U = cell(N,1);
for i = 1:N
    U{i} = randn(s(i),r);
end

load('init_U.mat')

[U,error_all,error_swa] = ADAM_swa(0,U,X);
% [U_full,error_full] = ADAM(2,U,X);

normX = norm(X);
maxiter = 20000;
figure
semilogy(error_all/normX)
%     semilogy(e1/normX)
hold on
semilogy(error_swa/normX)

legend('half precision','double precision')

xlabel('number of iterations')
ylabel('error')
title('d=[100,50,25], r=40')



