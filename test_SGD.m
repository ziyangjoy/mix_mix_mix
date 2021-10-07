%load aminoacids.mat
clear all;

addpath('./chop');
addpath('./tensor_toolbox-v3.2');

rng(12);

N = 3;


size_t = 80;
r = 30;

s = [size_t,size_t*2,size_t/2,size_t,size_t];


A = cell(N,1);
for i = 1:N
    A{i} = randn(s(i),r);
%     A{i} = 2*rand(s(i),r)-1;
%     A{i} = A{i}/max(A{i}(:));
end
X = ktensor(A);
X = double(tensor(X));


U = cell(N,1);
for i = 1:N
    U{i} = randn(s(i),r);
%     U{i} = U{i}/max(U{i}(:));
%     U{i} = ones(s(i),r);
end

% load('SGD_result_size50r20_double_newton.mat')
% U = U_result{88};

% [U_half,error_half] = SGD(2,U,X);
[U_full,error_full] = SGD_epoch(2,U,X);
% [U_half,error_half] = SGD_epoch(0,U,X);
% [U_half,error_half] = SGD_newsample_epoch(2,U,X);


% normX = norm(X);
% maxiter = 20000;
% figure
% semilogy(error_half(1:maxiter)/normX)
% 
% hold on
% semilogy(error_full(1:maxiter)/normX)
% 
% legend('half precision','double precision')
% 
% xlabel('number of iterations')
% ylabel('error')
% title('d=[100,50,25], r=40')



