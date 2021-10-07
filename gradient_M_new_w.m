function G = gradient_M_new_w(prec,U,X,ind)



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

num = cell(N,1);
M_num = cell(N,1);
for j = 1:N
   num{j} = zeros(s(j),1);
   [GC,GR] = groupcounts(ind(:,j));
   num{j}(GR) = GC;
   M_num{j}= ceil(M/s(j));
end

for i = 1:M
   ind_i = num2cell(ind(i,:));
   mat = [];
   for j = 1:N
       mat = [mat;U{j}(ind_i{j},:)] ;
   end
   v = sum(prod(mat))-X(ind_i{:});

   for j = 1:N
%        G{j}(ind_i{j},:) = G{j}(ind_i{j},:) + 2*prod(mat([1:j-1,j+1:end],:))*v/M_num{j};
%        G{j}(ind_i{j},:) = G{j}(ind_i{j},:) + 2*prod(mat([1:j-1,j+1:end],:))*v/(num{j}(ind_i{j})*s(j));
       G{j}(ind_i{j},:) = G{j}(ind_i{j},:) + 2*prod(mat([1:j-1,j+1:end],:))*v/num{j}(ind_i{j});
   end
end