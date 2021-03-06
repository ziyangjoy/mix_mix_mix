function x = gmres3(A,b,precf,precw,precr,iter_max,gtol)
dA = double(A);
sA = single(A);
hA = half(A);
dB = double(b);
sB = single(b);
hB = half(b);

if precf ~=0 && precf ~=1 && precf ~= 2, error('0 for half, 1 for single, 2 for double'), end
if precw ~=0 && precw ~=1 && precw ~= 2, error('0 for half, 1 for single, 2 for double'), end
if precr ~=1 && precr ~= 2, error('1 for single, 2 for double, dont use half for residue'), end

n = length(A);
realx = A\b;

if precf == 1
    %fprintf('**** uf = single.\n')
    ufs = 'single';
elseif precf == 2
    %fprintf('**** uf = double.\n')
    ufs = 'double';
else
    %fprintf('**** uf = half.\n')
    ufs = 'half';
end

if precw == 0
    %fprintf('**** uw = half.\n')
    uws = 'half';
    A = hA;
    b = hB;
    u = eps(half(1));
elseif precw == 1
    %fprintf('**** uw = single.\n')
    uws = 'single';
    A = sA;
    b = sB;
    u = eps('single');
    
else
    %fprintf('**** uw = double.\n')
    uws = 'double';
    A = double(A);
    b = dB;
    u = eps('double');
end

if precr == 1
    %fprintf('**** ur = single.\n')
    urs = 'single';
else 
    %fprintf('**** ur = double.\n')
    urs = 'double';
end

  

%LU factorization

if precf == 1
    [L,U,P] = lu(sA);
    LL = single(P'*L);
    x =  U\(L\(P*sB) );
elseif precf == 2
    [L,U,P] = lu(dA);
    LL = double(P'*L);
    x =  U\(L\(P*dB) );
else
    [L,U,p] = lu3(hA);
    I = half(eye(n)); 
    P = I(p,:);
    LL = half(P'*L);
    x =  half(U\(L\(double(P*hB))) );
end

%condition numbers for A

At = double(U)\double(L)\double(P)*dA;
kinfA = cond(dA,'inf');
kinfAt = cond(At,'inf');
condAx = norm(abs(inv(dA))*abs(dA)*abs(realx),inf)/norm(realx,inf);
condA = norm(abs(inv(dA))*abs(dA),'inf');

% initial x
 

if precw == 0
    x = half(x);
elseif precw == 2
    x = double(x);
else
    x = single(x);
end


iter = 0; dx = 0; rd = 0;

% # iterations
gmresits = [];

% residual norm
gmreserr = [];
ferr = inf(iter_max,1);
mu = inf(iter_max,1);
nbe = inf(iter_max,1);
cbe = inf(iter_max,1);
while true
    % copied from code in paper 
    % Compute size of errors, quantities in bounds
    ferr(iter+1) = double(norm(double(x)-realx,'inf')/norm(realx,'inf'));
    mu(iter+1) = norm(dA*(double(x)-realx),'inf')/(norm(dA,'inf')*norm(double(x)-realx,'inf')); 
    res = dB - dA*double(x);
    nbe(iter+1) = double(norm(res,'inf')/(norm(dA,'inf')*norm(double(x),'inf')+ norm(dB,'inf')));
    temp = double( abs(res) ./ (abs(dA)*abs(double(x)) + abs(dB)) );
    temp(isnan(temp)) = 0; % Set 0/0 to 0.
    cbe(iter+1) = max(temp);
    
    iter = iter + 1;
    if max([ferr(iter) nbe(iter) cbe(iter)]) <= u || iter > iter_max, break, end
    
    %Compute residual
    if precr == 1
        rd = sB - sA*single(x);
    elseif precr == 2
        rd = dB - dA*double(x);
    end
    
    %normalize residue 
    norm_rd = norm(rd,inf);
    rd1 = rd/norm_rd;
    
    if precw == 0
        [d, err, its, ~] = halfgmeres( A, half(zeros(n,1)), half(rd1), LL, U, n, 1, gtol);
    elseif precw == 1
        [d, err, its, ~] = gmres_sd( A, single(zeros(n,1)), single(rd1), LL, U, n, 1, gtol);
    else
        [d, err, its, ~] = gmres_sd( A, zeros(n,1), double(rd1), LL, U, n, 1, gtol);
    end
    
    % copied from code in paper 
    %Compute quantities in bounds for plotting
   
    lim(iter) = double( 2*u*cond(dA,'inf')*mu(iter));
    lim2(iter) = double(2*u*condA);  
    dact = dA\double(rd1);
    etai(iter) = norm(double(double(d)-dact),'inf')/norm(dact,'inf');   
    phi(iter) = min(lim(iter),lim2(iter))+etai(iter);
    
    %Record number of iterations gmres took
    gmresits = [gmresits,its];
    
    %Record final relative (preconditioned) residual norm in GMRES
    gmreserr = [gmreserr,err(end)];
    
    %Record relative (preconditioned) residual norm in each iteration of
    %GMRES (so we can look at convergence trajectories if need be)
    gmreserrvec{iter} = err;
    
    xold = x;
    
    %Update solution
    if precw == 0
        x = x + half(norm_rd)*half(d);
    elseif precw == 1
        x = x + single(norm_rd)*single(d);
    else
        x = x + norm_rd*double(d);
        
    end
    dx = norm(x-xold,'inf')/norm(x,'inf');
    
    %Check if dx contains infs, nans, or is 0
    if dx == Inf || isnan(double(dx))
        plt = 0;
        break;
    end
end