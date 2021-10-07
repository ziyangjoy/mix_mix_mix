function G = gradient_sto_new_w(prec,U,X,ind)



N = length(U);
s = size(X);

G = cell(N,1);
for j = 1:N
   G{j} = zeros(size(U{j}));
end

r = size(U{1},2);

fp.format = 'h';
fp.round = 5;

if prec == 0
    G = cellfun(@(x)chop(x,fp),G,'UniformOutput',0);
    U = cellfun(@(x)chop(x,fp),U,'UniformOutput',0);
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
   if prec == 0
       v = chop(v,fp);
   end

   for j = 1:N
       tmp = 2*prod(mat([1:j-1,j+1:end],:))*v/num{j}(ind_i{j});
       if prec == 0
           tmp = chop(tmp,fp);
           G{j}(ind_i{j},:) = chop(G{j}(ind_i{j},:) + tmp,fp);
       else
           G{j}(ind_i{j},:) = G{j}(ind_i{j},:) + tmp;
       end
   end
end