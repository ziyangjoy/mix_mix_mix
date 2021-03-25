function [x, error, its, flag] = halfgmeres( A, x, b, L, U, restrt, max_it, tol)

flag = 0;
its = 0;

%Ensure half working precision
A = half(A);
b = half(b);
x = half(x);
L = half(L);
U = half(U);


rtmp = b-A*x;

r = single(L)\single(rtmp);
r = single(U)\single(r);
r = half(r);

bnrm2 = norm(r );
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
error=[];

error(1) = norm( r ) / bnrm2;
if ( error(1) < tol ) return, end

[n,~] = size(A);                                  % initialize workspace
m = restrt;
V = half(zeros(n,m+1));
H = half(zeros(m+1,m));
cs = half(zeros(m,1));
sn = half(zeros(m,1));
e1    = half(zeros(n,1));
e1(1) = half(1.0);

for iter = 1:max_it,                              % begin iteration
    rtmp = half(b-A*x);
    r = single(U)\(single(L)\single(rtmp));
    r = half(r);
    
    V(:,1) = r / norm( r );
    s = norm( r )*e1;
    for i = 1:m,                     % construct orthonormal basis via GS
        its = its+1;
        vcur = V(:,i);      
        
        vcur = single(U)\(single(L)\(single(A)*single(vcur)));
        
        w = half(vcur);
      
        for k = 1:i,
            H(k,i)= w'*V(:,k);
            w = w - H(k,i)*V(:,k);
        end
        H(i+1,i) = norm( w );
        V(:,i+1) = w / H(i+1,i);
        for k = 1:i-1,                              % apply Givens rotation
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
        end
        [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
        temp   = cs(i)*s(i);                        % approximate residual norm
        s(i+1) = -sn(i)*s(i);
        s(i)   = temp;
        H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
        error((iter-1)*m+i+1)  = abs(s(i+1)) / bnrm2;
        if ( error((iter-1)*m+i+1) <= tol ),                        % update approximation
            y = H(1:i,1:i) \ s(1:i);                 % and exit
            addvec = V(:,1:i)*y;
            x = x + addvec;
            break;
        end
    end
    
    if ( error(end) <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    addvec = V(:,1:m)*y;
    x = x + addvec;                            % update approximation
    rtmp = b-A*x;
    r = single(U)\(single(L)\single(rtmp));           % compute residual
    r = half(r);
    s(i+1) = norm(r);
    error = [error,s(i+1) / bnrm2];                        % check convergence
    if ( error(end) <= tol ), break, end;
end

if ( error(end) > tol ) flag = 1; end;                 % converged




function [ c, s ] = rotmat( a, b )
%
% Compute the Givens rotation matrix parameters for a and b.
%
if ( b == 0.0 ),
    c = 1.0;
    s = 0.0;
elseif ( abs(b) > abs(a) ),
    temp = a / b;
    temp2 = temp*temp; 
    s = 1.0 / half(sqrt(single( 1.0 + temp2 )));
    c = temp * s;
else
    temp = b / a;
    temp2 = temp*temp;
    c = 1.0 / half(sqrt(single( 1.0 + temp2 )));
    s = temp * c;
end