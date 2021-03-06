%cp_als_new(X,R,method,uf,uw,ur,uk,giters,gtol,varargin)
% parameters: 
% X: tensor; R: rand;
% method: 0 for least square by gmres; 1 for the original method
% uf, uw, ur, uk: precision for factorization, working, residue, khatri-rao
% giters: max number of iterations for gmres
% gtol: tolerance for gmres
% Other variables:      
%      'tol' - Tolerance on difference in fit {1.0e-4}
%      'maxiters' - Maximum number of iterations {50}
%      'dimorder' - Order to loop through dimensions {1:ndims(A)}
%      'init' - Initial guess [{'random'}|'nvecs'|cell array]
%      'printitn' - Print fit every n iterations; 0 for no printing {1}
%      'fixsigns' - Call fixsigns at end of iterations {true}

% load the tensor for test
load aminoacids.mat
% test 1: using different precision for khatri-rao product: 
% error is measured by 1-fit.
max_iters = [1,5,10,20,50];
error_list_0 = [];
error_list_2 = [];
for i = max_iters
    [~,~,~,fit] = cp_als_new(X,50,0,2,2,2,0,3,0.0001,'maxiters',i,'printitn',0);
    error0 = 1-fit;
    error_list_0(end+1)=error0;
    [~,~,~,fit2] = cp_als_new(X,50,0,2,2,2,2,3,0.0001,'maxiters',i,'printitn',0);
    error2 = 1-fit2;
    error_list_2(end+1)=error2;
end
% number of max iterations versus error
figure
title('different khatri-rao product precision')
plot(max_iters,error_list_0,'color','r')
hold on
plot(max_iters,error_list_2,'color','g')

% from test 1 we can see that using different precision for khatri-rao
% product doesn't affect the convergence a lot.

% test 2: using different precision for factorization

% for i = max_iters
%     [~,~,~,fit] = cp_als_new(X,20,0,0,2,2,2,3,0.0001,'maxiters',i,'printitn',0);
%     error0 = 1-fit;
%     error_list_0(end+1)=error0;
%     [~,~,~,fit2] = cp_als_new(X,20,0,2,2,2,2,3,0.0001,'maxiters',i,'printitn',0);
%     error2 = 1-fit2;
%     error_list_2(end+1)=error2;
% end
% 
% % number of max iterations versus error
% % low factorization precision sometimes leads to bad gmres result
% 
% figure
% title('different factorization precision')
% plot(max_iters,error_list_0,'color','r')
% hold on
% plot(max_iters,error_list_2,'color','g')



