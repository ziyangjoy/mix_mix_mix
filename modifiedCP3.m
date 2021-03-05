function [P,Uinit,output,error] = modifiedCP3(X,R,method,uf,uw,ur,uk,max_iter_gmres,varargin)
%CP_ALS Compute a CP decomposition of any type of tensor.
%
%   P = CP_ALS(X,R) computes an estimate of the best rank-R
%   CP model of a tensor X using an alternating least-squares
%   algorithm.  The input X can be a tensor, sptensor, ktensor, or
%   ttensor. The result P is a ktensor.
%
%   P = CP_ALS(X,R,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'tol' - Tolerance on difference in fit {1.0e-4}
%      'maxiters' - Maximum number of iterations {50}
%      'dimorder' - Order to loop through dimensions {1:ndims(A)}
%      'init' - Initial guess [{'random'}|'nvecs'|cell array]
%      'printitn' - Print fit every n iterations {1}
%
%   [P,U0] = CP_ALS(...) also returns the initial guess.
%
%   [P,U0,out] = CP_ALS(...) also returns additional output that contains
%   the input parameters.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   P = cp_als(X,2);
%   P = cp_als(X,2,'dimorder',[3 2 1]);
%   P = cp_als(X,2,'dimorder',[3 2 1],'init','nvecs');
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
%   [P,U0,out] = cp_als(X,2,'dimorder',[3 2 1],'init',U0);
%   P = cp_als(X,2,out.params); %<-- Same params as previous run
%
%   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda.
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: cp_als.m,v 1.7 2010/03/19 23:46:32 tgkolda Exp $


%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParamValue('tol',1e-4,@isscalar);
params.addParamValue('maxiters',50,@(x) isscalar(x) & x > 0);
params.addParamValue('dimorder',1:N,@(x) isequal(sort(x),1:N));
params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addParamValue('printitn',1,@isscalar);
params.parse(varargin{:});

%% Copy from params object
fitchangetol = params.Results.tol;
maxiters = params.Results.maxiters;
dimorder = params.Results.dimorder;
init = params.Results.init;
printitn = params.Results.printitn;

%% Error checking 

%% Set up and error checking on initial guess for U.
if iscell(init)
    Uinit = init;
    if numel(Uinit) ~= N
        error('OPTS.init does not have %d cells',N);
    end
    for n = dimorder(2:end);
        if ~isequal(size(Uinit{n}),[size(X,n) R])
            error('OPTS.init{%d} is the wrong size',n);
        end
    end
else
    % Observe that we don't need to calculate an initial guess for the
    % first index in dimorder because that will be solved for in the first
    % inner iteration.
    if strcmp(init,'random')
        Uinit = cell(N,1);
        for n = dimorder(1:end)
            Uinit{n} = rand(size(X,n),R);
        end
    elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
        Uinit = cell(N,1);
        for n = dimorder(1:end)
            Uinit{n} = nvecs(X,n,R);
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;
if printitn>0
  fprintf('\nCP_ALS:\n');
end

%% Main Loop: Iterate until convergence
Xnlist = {};
Xcontent = X.data;
for n = 1:N
    Xnlist{end+1}=n_mode_unfold(Xcontent,n);
end

for iter = 1:maxiters   
    fitold = fit;
    for n = dimorder(1:end)
        if method == 1
            Unew = mttkrp(X,U,n);
            Y = ones(R,R);
            for i = [1:n-1,n+1:N]
                Y = Y .* (U{i}'*U{i});
            end
            Unew = (Y \ Unew')';
        elseif method == 0 
            Kr_product = kr3_all_except_mode(U,n,uk);
            Xn = Xnlist{n};
            Xnt = Xn';
            Ant = getx(Kr_product,Xnt,uf,uw,ur,max_iter_gmres,0.0001);
            Unew = Ant';
        end
    end
            
    % Normalize each vector to prevent singularities in coefmatrix
    if iter == 1
        lambda = sqrt(sum(Unew.^2,1))'; %2-norm
    else
        lambda = max( max(Unew,[],1), 1 )'; %max-norm
    end
    Unew = Unew * spdiags(1./lambda,0,R,R);
    if issparse(Unew)
        U{n} = full(Unew);   % for the case R=1
    else
        U{n} = Unew;
    end

    P = ktensor(lambda,U);
    normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
    fit = 1 - (normresidual / normX); %fraction explained by model
    fitchange = abs(fitold - fit);

    if mod(iter,printitn)==0
      fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
    end

    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end

end

%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);
% Fix the signs
P = fixsigns(P);
kproduct = kr3_all_except_mode(P,0,uk);
%error = norm(kproduct-Xnlist{end});
if printitn>0
  normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
  fit = 1 - (normresidual / normX); %fraction explained by model
  fprintf(' Final fit = %e \n', fit);
end
error = sqrt(normX-norm(P));
output = struct;
output.params = params.Results;