function X = norm3(A,p)
if nargin < 2 || p == f
    flatA = reshape(A,[],1);
    AA = flatA.*flatA;
    X = sqrt(sum(AA));
    
end
     
