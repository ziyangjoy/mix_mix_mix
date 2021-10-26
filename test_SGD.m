%load aminoacids.mat
clear all;

addpath('./chop');
addpath('./tensor_toolbox-v3.2');

rng(41);

N = 3;


size_t = 50;
r = 30;

% s = [80,50,30,20];
s = [size_t,size_t,size_t,size_t,size_t];

fp.format = 'c';
fp.params = [4,7] ;

fp.round = 5;

A = cell(N,1);
for i = 1:N
    A{i} = randn(s(i),r);
    A{i} = chop(A{i},fp);
%     A{i} = A{i}./vecnorm(A{i});
%       A{i} = (randi(7,s(i),r) - 4)/3;
%     A{i} = 2*rand(s(i),r)-1;
%     A{i} = A{i}/max(A{i}(:));
end
X = ktensor(A);
X = double(tensor(X));
X = X/max(abs(X(:)))*15;

X = X + 0.000*randn(s(1:N));


U = cell(N,1);

% w = 1;
% r = r + 5;
for i = 1:N
%     U{i} = A{i} + 0.05*randn(s(i),r);
    U{i} = randn(s(i),r)/2;
%     w = w.*vecnorm(U{i}).'/5;
%     U{i} = U{i}./vecnorm(U{i})*5;
%     U{i} = U{i}/max(U{i}(:));
%     U{i} = ones(s(i),r);
end


% load('SGD_result_size50r20_double_newton.mat')
% U = U_result{88};

% [U_half,error_half] = SGD(2,U,X);
% [U_full,error_full] = SGD_epoch_sto(0,U,X);
[U_full,error_full] = ADAM_epoch_unbiased(0,U,X);
% [U_full,error_full] = ADAM_epoch(2,U,X);

% [U_full,error_full] = SGD_epoch(2,U,X);
% [U_full,error_full] = SVRG_epoch(2,U,X);
% [U_half,error_half] = SGD_epoch_weight_decay(2,U,X);
% [U_half,error_half] = SGD_newsample_epoch(2,U,X);


normX = norm(X(:));
maxiter = 20000;
figure
semilogy(error_half/normX)

hold on
semilogy(error_full/normX)

legend('half precision','double precision')

xlabel('number of epoch')
ylabel('error')
title('d=[60,60,60], r=40')



