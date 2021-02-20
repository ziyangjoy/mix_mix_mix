error_norm = [];
error_product = [];
for size = [10,20,30,40,50]
    A = size*rand(size);
    hA = half(A);
    en = norm(double(A-hA),'fro');
    invhA = inv(double(half(A)));
    ep = norm(A*invhA-eye(size),'fro');
    error_norm(end+1)=en;
    error_product(end+1)=ep;
end
plot(error_norm,error_product)
    