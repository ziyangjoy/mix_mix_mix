function G = gradient_full_ownfunc(prec,U,X)

N = length(U);

G = cell(N,1);
for j = 1:N
   G{j} = (zeros(size(U{j})));
end

r = size(U{1},2);

fp.format = 'h';
fp.round = 5;

if prec == 0
    U = cellfun(@(x)half(x),U,'UniformOutput',0);
    G = cellfun(@(x)half(x),G,'UniformOutput',0);
end

T = gen_ten(U);


for j = 1:N
    M_X = tenmode_k(X,j);
    M_T = tenmode_k(T,j);
    tmp = {U{[N:-1:j+1,j-1:-1:1]}};
    V = khatrirao_Z(tmp);
    
    P = M_T - M_X;
%     G{j} = P.'*V*inv(V'*V);
%     G{j} = P.'*V*diag(1./vecnorm(V).^2);
    G{j} = 2*P.'*V;
    
    
end