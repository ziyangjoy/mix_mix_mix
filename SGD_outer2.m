function [u,e,T] = SGD_outer2(prec,initial, X,iterCG,iterSG,tol,R,step)
s = size(X);
N = ndims(X);
% A = cell(N,1);
A = initial;
e = [];
% for i = 1:N
%      A{i} = randn(s(i),R);
% end
tic;
T = zeros(1,iterSG);
for i = 1:iterSG
    [A,error] = CG_cp2(prec,X,R,A,iterCG,tol,step);
    e(end+1)=error;
    T(i)=toc;
end
u = A;
normX = norm(X);
% e = e/normX;