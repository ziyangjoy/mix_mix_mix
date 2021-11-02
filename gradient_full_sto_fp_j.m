function G = gradient_full_sto_fp_j(prec,U,X,fp,j)

N = length(U);
s = cellfun(@(x)size(x,1),U).';



r = size(U{1},2);

fp_sum.format = 'h';


% fp.format = 'h';
% fp.format = 'c';
% fp.params = [4,7] ;
% fp.params = [5,3] ;

% fp.round = 5;
U_old = U;
if prec == 0
%     G = cellfun(@(x)chop(x,fp),G,'UniformOutput',0);
    U = cellfun(@(x)chop(x,fp),U,'UniformOutput',0);
%     U = cellfun(@(x)half(x),U,'UniformOutput',0);
%     G = cellfun(@(x)half(x),G,'UniformOutput',0);
end

% T = ktensor(U);
T = double(tensor(ktensor(U)));
if prec == 0
    V = ones(1,r);
    for k = N:-1:1
        V = khatrirao(V,U{k});
        if prec == 0
           V = chop(V); 
        end
    end
    T = chop(reshape(sum(V,2),s),fp);
end


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
   G = M_T - M_X;
%    G = chop((M_T-M_X),fp);
   G = chop(2*chop(G.'*V,fp),fp);

%        G{j} = chop(chop(2*M_T.'*V,fp) - chop(2*M_X.'*V,fp),fp);
else
   G = 2*(M_T-M_X).'*V;
end
