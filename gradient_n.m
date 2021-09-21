function G = gradient_n(prec,U,X,n)

N = length(U);

G = cell(N,1);
for j = 1:N
   G{j} = (zeros(size(U{j})));
end

r = size(U{1},2);

if prec == 0
    U = cellfun(@(x)half(x),U,'UniformOutput',0);
    G = cellfun(@(x)half(x),G,'UniformOutput',0);
end

ind_i = num2cell(n);
mat = [];
for j = 1:N
   mat = [mat;U{j}(ind_i{j},:)] ;
end
p = prod(mat);
v = sum(p)-X(ind_i{:});
for j = 1:N
   G{j}(ind_i{j},:) = 2*p./U{j}(ind_i{j},:)*v;
end
