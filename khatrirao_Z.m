function M = khatrirao_Z(U)

N = length(U);
r = size(U{1},2);
D = prod(cellfun(@(x)size(x,1),U));
M = zeros(D,r);
% M = [];
for i = 1:r
    tmp = 1;
    for j = 1:N
       tmp = kron(tmp,U{j}(:,i)) ;
    end
    M(:,i) = tmp;
end


