R = 20;
N = ndims(X);
dimorder = [1:N];
normX = norm(X);
Uinit = cell(N,1);
for n = dimorder(1:end)
    Uinit{n} = rand(size(X,n),R);
end
U = Uinit;
Xnlist = {};
for n = dimorder(1:end)
    Xnlist{end+1}=n_mode_unfold(Xcontent,n);
end
Xcontent = X.data;
maxiters = 1;
uk = 2;
uf = 0;
uw = 2;
ur = 2;
giters = 1;
gtol = 0.0001;
printitn = 0;
fit = 0;
fitchangetol = 0.00001;
for iter = 1:maxiters
    fitold = fit;
    for n = dimorder(1)
        Kr_product = kr3_for_cp(U,n,uk);
        Xn = Xnlist{n};
        Xnt = Xn';
        Ant = getx(Kr_product,Xnt,uf,uw,ur,giters,gtol);
        Unew = Ant';
        if n == dimorder(end)
            U_mttkrp = Unew;
        end
        if iter == 1
            lambda = sqrt(sum(Unew.^2,1))'; %2-norm
        else
            lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
        end
        
        Unew = bsxfun(@rdivide, Unew, lambda');
        
        U{n} = Unew;
        UtU(:,:,n) = U{n}'*U{n};
    end
    
    P = ktensor(lambda,U);
    
    % This is equivalent to innerprod(X,P).
    
    iprod = sum(sum(P.U{dimorder(end)} .* U_mttkrp) .* lambda');
    if normX == 0
        fit = norm(P)^2 - 2 * iprod;
    else
        normresidual = sqrt( normX^2 + norm(P)^2 - 2 * iprod );
        fit = 1 - (normresidual / normX); %fraction explained by model
    end
    fitchange = abs(fitold - fit);
    
    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        flag = 0;
    else
        flag = 1;
    end
    
    if (mod(iter,printitn)==0) || ((printitn>0) && (flag==0))
        fprintf(' Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange);
    end
    
    % Check for convergence
    if (flag == 0)
        break;
    end
end