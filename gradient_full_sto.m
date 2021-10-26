function G = gradient_full_sto(prec,U,X)

N = length(U);

G = cell(N,1);
for j = 1:N
   G{j} = (zeros(size(U{j})));
end

r = size(U{1},2);

fp.format = 'h';
fp.format = 'c';
fp.params = [4,7] ;
% fp.params = [5,3] ;

fp.round = 5;

if prec == 0
    G = cellfun(@(x)chop(x,fp),G,'UniformOutput',0);
    U = cellfun(@(x)chop(x,fp),U,'UniformOutput',0);
%     U = cellfun(@(x)half(x),U,'UniformOutput',0);
%     G = cellfun(@(x)half(x),G,'UniformOutput',0);
end

% T = ktensor(U);
T = double(tensor(ktensor(U)));

for j = 1:N
    M_X = tenmode_k(X,j);
    M_T = tenmode_k(T,j);
%     M_X = double(tenmat(X,[1:j-1,j+1:N],[j]));
%     M_T = double(tenmat(T,[1:j-1,j+1:N],[j]));
    V = ones(1,r);
    for k = [N:-1:j+1,j-1:-1:1]
        V = khatrirao(V,U{k});
        if prec == 0
           V = chop(V); 
        end
    end
    if prec == 0
       M_X = chop(M_X,fp); 
       M_T = chop(M_T,fp); 
       G{j} = chop(2*(M_T-M_X),fp);
       G{j} = chop(G{j}.'*V,fp);
%        G{j} = chop(2*(M_T-M_X).'*V,fp);
    else
       G{j} = 2*(M_T-M_X).'*V;
    end
    
end