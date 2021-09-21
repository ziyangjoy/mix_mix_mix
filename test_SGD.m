%load aminoacids.mat
N = 3;

size = 10;

s = [size/2,size*2,size,size/2,size/2];
r = 10;
A = cell(N,1);
for i = 1:N
    A{i} = randn(s(i),r);
end
X = ktensor(A);


U = cell(N,1);

for i = 1:N
    U{i} = randn(s(i),r);
end

prec = 1;
M = 100;
alpha = 0.01;

for i = 1:100
    G = gradient_M(prec,U,X,M);
    for j = 1:N
        U{j} = U{j} - alpha*G{j};
    end
end

U = cellfun(@(x)double(x),U,'UniformOutput',0);
nX = ktensor(U);

error = norm(minus(full(nX),full(X)));