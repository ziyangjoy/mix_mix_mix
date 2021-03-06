function x = getx(A,B,uf,uw,ur,max_iter,gtol)
[nm,J] = size(B);
[mn,R] = size(A);
if nm ~= mn
    error(message('input matrix dimension dont match'));
end
x = zeros(R,J);
if mn == R
    for i = 1:J
        bi = B(:,i);
        xi = gmres3(A,bi,uf,uw,ur,max_iter,gtol);
        x(:,i)=xi;
    end
else
    At = A';
    AtA = At*A;
    AtB = At*B;
    for i = 1:J
        bii = AtB(:,i);
        xi = gmres3(AtA,bii,uf,uw,ur,max_iter,gtol);
        x(:,i)=xi;
    end
end