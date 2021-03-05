% error is defined by sqrt(norm(P)-norm(X)), where P is the final output. 
% modifiedCP3(X,r,method,uf,uw,ur,uk,max_iter_gmres,varargin)
% X: tensor
% r: rank
% method: 0 represent gmres methods written my Jiajia Wang; 1 represent the
% original cp-als method in tensor toolbox
% uf: factorization precision; uw: working precision; ur: residual
% precision; uk: khatri rao product precision
% max_iter_gmres: max number of iterations for gmres
% varargin: 
%      'tol' - Tolerance on difference in fit {1.0e-4}
%      'maxiters' - Maximum number of iterations {50}
%      'dimorder' - Order to loop through dimensions {1:ndims(A)}
%      'init' - Initial guess [{'random'}|'nvecs'|cell array]
%      'printitn' - Print fit every n iterations {1}

X = rand(200,100,25,50);
X = tensor(X);
errorlist_0 = [];
errorlist_2 = [];
rank = [10,20,30,40,50,200];
% by change the maxiters, we can see the speed of convergence of cp-als
% with different precision of khatri-rao product.

for r = [10,20,30,40,50,200]
    [~,~,~,error0] = modifiedCP3(X,r,0,0,2,2,0,2,'maxiters',1);
    errorlist_0(end+1)=error0;
    [~,~,~,error2] = modifiedCP3(X,r,0,0,2,2,2,2,'maxiters',1);
    errorlist_2(end+1)=error2;
end
figure
plot(rank,errorlist_0)
hold
plot(rank,errorlist_2)
