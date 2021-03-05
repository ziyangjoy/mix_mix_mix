function P = kr3(A,B,prec)
[ma,na] = size(A);
[mb,nb] = size(B);
if prec == 0
    A = half(A);
    B = half(B);
elseif prec == 1
    A = single(A);
    B = single(B);
else
    A = double(A);
    B = double(B);
end
if na ~= nb
    error(message('input matrix dimension dont match'));
end
P = zeros(ma*mb,na);
for i = (1:na)
    P(:,i) = kron(A(:,i),B(:,i));
end
