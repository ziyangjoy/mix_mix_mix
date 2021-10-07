%load aminoacids.mat
clear all;

addpath('./chop');
addpath('./tensor_toolbox-v3.2');

rng(15);

N = 3;

prec = 0;

size = 50;
r = 20;

s = [size/2,size*2,size];

A = cell(N,1);
for i = 1:N
    A{i} = double(half(randn(s(i),r)));
end
X = ktensor(A);


U = cell(N,1);
for i = 1:N
    U{i} = randn(s(i),r)+.1;
end

% [U_half,error_half] = ADAM(0,U,X);
[U_full,error_full] = ADAM(2,U,X);

normX = norm(X);
maxiter = 20000;
figure
semilogy(error_half(1:maxiter)/normX)
%     semilogy(e1/normX)
hold on
semilogy(error_full(1:maxiter)/normX)

legend('half precision','double precision')

xlabel('number of iterations')
ylabel('error')
title('d=[100,50,25], r=40')



