function V = mymttkrp1(prec,X,U,n)

N = ndims(X);

if n == 1
    R = size(U{2},2);
else
    R = size(U{1},2);
end
   
%% Computation
szl = prod(size(X,1:n-1)); 
szr = prod(size(X,n+1:N));
szn = size(X,n);
if prec == 2
    if n == 1
        Ur = mykhatrirao(2,U{2:N},'r');
        Y = reshape(X.data,szn,szr);
        V =  Y * Ur;
    elseif n == N
        Ul = mykhatrirao(2, U{1:N-1},'r');
        Y = reshape(X.data,szl,szn);
        V = Y' * Ul;
    else
        Ul = mykhatrirao(2, U{n+1:N},'r');
        Ur = reshape(mykhatrirao(2, U{1:n-1},'r'), szl, 1, R);
        Y = reshape(X.data,[],szr);
        Y = Y * Ul;
        Y = reshape(Y,szl,szn,R);
        V = zeros(szn,R);
        for r =1:R
            V(:,r) = Y(:,:,r)'*Ur(:,:,r);
        end
    end
elseif prec == 1
    if n == 1
        Ur = mykhatrirao(1, U{2:N},'r');
        Y = single(reshape(X.data,szn,szr));
        V =  Y * Ur;
    elseif n == N
        Ul = mykhatrirao(1, U{1:N-1},'r');
        Y = single(reshape(X.data,szl,szn));
        V = Y' * Ul;
    else
        Ul = mykhatrirao(1, U{n+1:N},'r');
        Ur = reshape(mykhatrirao(1, U{1:n-1},'r'), szl, 1, R);
        Y = single(reshape(X.data,[],szr));
        Y = Y * Ul;
        Y = reshape(Y,szl,szn,R);
        V = single(zeros(szn,R));
        for r =1:R
            V(:,r) = Y(:,:,r)'*Ur(:,:,r);
        end
    end
elseif prec == 0
    if n == 1
        Ur = mykhatrirao(0, U{2:N},'r');
        Y = half(reshape(X.data,szn,szr));
        V =  Y * Ur;
    elseif n == N
        Ul = mykhatrirao(0, U{1:N-1},'r');
        Y = half(reshape(X.data,szl,szn));
        V = Y' * Ul;
    else
        Ul = mykhatrirao(0, U{n+1:N},'r');
        Ur = reshape(mykhatrirao(0, U{1:n-1},'r'), szl, 1, R);
        Y = half(reshape(X.data,[],szr));
        Y = Y * Ul;
        Y = reshape(Y,szl,szn,R);
        V = half(zeros(szn,R));
        for r =1:R
            V(:,r) = Y(:,:,r)'*Ur(:,:,r);
        end
    end
    
end


