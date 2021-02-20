function P = kr3(A,B)
[ma,na] = size(A);
[mb,nb] = size(B);
if na ~= nb
    error(message('input matrix dimension dont match'));
end
P = zeros(ma*mb,na);
for i = (1:na)
    P(:,i) = kron(A(:,i),B(:,i));
end
