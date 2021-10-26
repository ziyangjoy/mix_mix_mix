function G = gradient_full_new(prec,U,X)

N = length(U);

G = cell(N,1);
for j = 1:N
   G{j} = (zeros(size(U{j})));
end

r = size(U{1},2);


if prec == 0
    U = cellfun(@(x)half(x),U,'UniformOutput',0);
%     U = cellfun(@(x)double(half(x)),U,'UniformOutput',0);
    G = cellfun(@(x)half(x),G,'UniformOutput',0);
    T = gen_ten(U);
else
    T = double(tensor(ktensor(U)));
end


% T = gen_ten(U);
% T = double(tensor(ktensor(U)));
% T = ktensor(U);

for j = 1:N
    M_X = tenmode_k(X,j);
    M_T = tenmode_k(T,j);
%     M_X = double(tenmat(X,[1:j-1,j+1:N],[j]));
%     M_T = double(tenmat(T,[1:j-1,j+1:N],[j]));
    
    tmp = U([N:-1:j+1,j-1:-1:1]);
%     V = khatrirao_Z(tmp);
    if prec == 0
        V = khatrirao_Z(tmp);
    else
        V = ones(1,r);
        for k = [N:-1:j+1,j-1:-1:1]
            V = khatrirao(V,U{k});
        end
    end
%     
%     if prec == 0
%        V = half(V);
% %        M_X = half(M_X);
%        M_T = half(M_T);
%     end
    
    P = M_T - M_X;
%     G{j} = P.'*V*inv(V'*V);
%     G{j} = P.'*V*diag(1./vecnorm(V).^2);
    G{j} = 2*P.'*V;
    
%     if prec == 0
%        G{j} = half(G{j}); 
%     end
    
end