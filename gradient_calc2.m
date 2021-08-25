function G = gradient_calc2(prec,X,U)
s = prod(size(X));
%s = 1;
N = ndims(X);
U = cellfun(@(x)(double(x)),U,'UniformOutput',0);
Uk = ktensor(U);
Y = minus(full(X),full(Uk));
Y = Y./s;
G = cell(N,1);
for i = 1:N
    G{i} = mymttkrp1(prec,Y,U,i);
end
G = cellfun(@(x)(-x),G,'UniformOutput',0);
