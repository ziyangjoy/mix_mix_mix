function V = gradient_n(prec,U,X,M)

N = length(U);
s = size(X);

r = size(U{1},2);

if prec == 0
    U = cellfun(@(x)half(x),U,'UniformOutput',0);
end

ind = [];
for i = 1:N
    tmp = randi(s(i),M,1);
    ind = [ind,tmp];
end
mat = [];
for i = 1:M
   ind_i = mat2cell(ind(i,:));
   mat = [];
   for j = 1:N
       mat = [mat;U{j}(ind_i{j},:)] ;
   end
   v = sum(prod(mat))-X(ind_i{:});
   
end