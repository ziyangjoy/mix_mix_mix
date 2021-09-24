function G = gradient_M(prec,U,X,ind)

N = length(U);
s = size(X);

G = cell(N,1);
for j = 1:N
   G{j} = zeros(size(U{j}));
end

r = size(U{1},2);

if prec == 0
    G = cellfun(@(x)half(x),G,'UniformOutput',0);
    U = cellfun(@(x)half(x),U,'UniformOutput',0);
end

M = size(ind,1);

for i = 1:M
   ind_i = num2cell(ind(i,:));
   mat = [];
   for j = 1:N
       if prec == -1
            mat = [mat;half(U{j}(ind_i{j},:))] ;
       else
           mat = [mat;U{j}(ind_i{j},:)] ;
       end
   end
   v = sum(prod(mat))-X(ind_i{:});
%    p = prod(mat);

   for j = 1:N
       G{j}(ind_i{j},:) = G{j}(ind_i{j},:) + 2*prod(mat([1:j-1,j+1:end],:))*v/M;
   end
end