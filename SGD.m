function [U,error_all] = SGD(prec,U,X)

% rng(15);

s = size(X);
M = floor(prod(s)/500);
M = 300;


alpha = .01;

N = length(s);
r = size(U{1},2);



error_all = [];
v = cell(N,1);
eta = 0.9;

for j = 1:N
   v{j} = zeros(size(U{j})) ;
end

for t = 1:50000
    n = [];
    for j = 1:N
        tmp = randi(s(j),M,1);
        n = [n,tmp];
    end
    G = gradient_M_new_w(prec,U,X,n);

    for j = 1:N
       v{j} = alpha*(G{j}) + eta*v{j};
       U{j} = U{j} - v{j};
    end
    nU = cellfun(@(x)double(x),U,'UniformOutput',0);
    nX = ktensor(nU);

    error = norm(minus(full(X),full(nX)));
    error_all = [error_all,error];
    
    if mod(t,500)==0
       error, 
    end
    

    
end