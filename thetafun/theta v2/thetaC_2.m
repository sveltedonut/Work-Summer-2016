function [theta_out, Hold, DHold, D2Hold, iter] = ...
                     thetaC(theta_in, Wfdobj, dataStr, Slambda, ...
                              itermax, crit)
%  This computes a score minimizing the regularized negative log
%  likelihood for a response pattern in U.  thata_in is the initial
%  value, and functional data object WFDOBJ is an expansion for the
%  functions log[P/(1-P)] associated with each item.  
%  The negative log likelihood for each theta value is regularized by 
%  attaching a penalty proportional to the squared difference
%  between the theta value and the icorresponding quantile of the scaled 
%  beta distribution over [0,n] with parameters ALPHA, BETA, H0 and HN.  
%  Only the quantile values for each case are required, and these are
%  computed externally and supplied in argument THETA_Q.
%  
%  Arguments are:
%  THETA_IN ... Initial values for abilities in [0,n].  It is assumed
%               that values corresponding to rows of U that are all
%               0's will be 0, and values corresponding to all 1's will be
%               n's.
%               This is a column vector of size N, where N is the number
%               of examinees.  The numbers will be in [0,n] and will be
%               floating point numbers rather than integers, where n is
%               the number of test items.
%  WFDOBJ   ... Functional data object defining the log-odds functions
%               for item characteristic curves and their smoothing.
%               A functional data object, consisting of n curves, each
%               a common basis object and a set of coefficients for
%               the expansion of the function in terms of the basis
%               functions.  The coefficient matrix (retrieved by the
%               getcoefs(basisobj) command) will have n columns and
%               K rows, where K is the number of basis functions.  In this
%               we are using spline basis functions.
%  DATASTR  ... A struct object with fields
%               U   ... N by n matrix of 0's and 1's indicating item 
%                       failure and success, respectively, for each of N  
%                       examinees and for each of n items.
%               IND ... The indices of the rows of U which are neither 
%                       all 0's or all 1's.  These numbers will be integers
%                       in the range 1 ... N and the size of the one-col
%                       array will be the numbere of examinees having 
%                       scores that are neither 0 or n.
%               TQ  ... A vector of length N containing the quantiles of 
%                       the tilted scaled beta distribution corresponding
%                       to the sum scores.  
%                       This is a single column of N float numbers.
%                       If SLAMBDA is missing or set to zero, these values  
%                       are not used and the argument can be [].
%  SLAMBDA  ... A smoothing parameter used to shrink ability parameters in
%               THETA towards the corresponding quantile of the 
%               scaled beta distribution fit to the sum scores.
%               This is a non-negative float number.
%  ITERMAX  ... Maximum number of iterations for computing the optimal
%               theta values.  Default is 20.
%               This is a nonnegative integer.
%  CRIT     ... Criterion for convergence of optimization.  Default is 1e-8
%               This is a positive float.
%
%  Last modified 3 January 2014 by Jim Ramsay

%  set dfault argument values

if nargin < 3
    error('Less than three arguments supplied.');
end

if nargin < 6,  crit    = 1e-8;  end 
if nargin < 5,  itermax = 20;    end 
if nargin < 4,  Slambda =  0;    end 

%  extract data objects from DATASTR 

U    = dataStr.U;
tQ   = dataStr.tQ;
ind  = dataStr.ind;
Nin  = sum(ind);
Uin  = U(ind,:);
tQin = tQ(ind);

[N,n] = size(U); 

tin = theta_in(ind);

if itermax == 0 
    theta_out = theta_in;
    return;  
end

if Slambda > 0
    tQ = tQ(:);
    if any(size(theta_in) ~= size(tQ))
        error('SLAMBDA is positive but size of THETA_Q is incorrect.');
    end
end

%  retrieve variables from fitstruct

wbasis  = getbasis(Wfdobj);    % basis for W functions

%  indices of scores on the boundaries and in the interior

scrrng  = getbasisrange(wbasis);
scrmin  = scrrng(1);
scrmax  = scrrng(2);

% compute initial probability values for interior scores

%% theta_2_0_0
[iter, tina, active, Ua, Hold] = theta_2_0_0(tin, Wfdobj, Uin, N, Slambda, tQin, crit);

while ~isempty(active) && iter < itermax
    %% theta_2_1_0
    
    [iter, tinanew, Hanew] = theta_2_1_0(iter, tina, Wfdobj, N, Ua, Slambda, tQin, active, scrmin, scrmax, tin, Hold, crit);
    
%     [iter, step, tinanew, Hanew, failindex] = theta_1_1_0(iter, tina, Wfdobj, N, Ua, Slambda, tQin, active, scrmin, scrmax, tin, Hold);
    
%     [iter, step] = theta_0_1_0(iter, tina, Wfdobj, N, Ua, Slambda, tQin, active);
%     
%     [tinanew, Hanew, fac, stepiter, failindex] = theta_0_1_1(tina, step, scrmin, scrmax, tin, active, Hold);
   
    
%     [tinanew, Hanew] = theta_1_2_1(step, scrmin, scrmax, tinanew, active, Hanew, failindex, Wfdobj, Ua, Slambda, tQin, N, Hold, crit);
%     while sum(failindex) > 0 && stepiter < 10
%         
%         
% %         
% %         [stepiter, tinanew] = theta_0_2_0(stepiter, step, scrmin, scrmax, failindex, fac, tinanew, active);
% %         
% %         [Hanew, failindex, fac] = theta_0_2_1(Hanew, failindex, fac, Wfdobj, tinanew, Ua, Slambda, tQin, N, Hold, crit, active);
%         
%     end
    %% theta_2_1_1
    
    [tina, tin, active, Hold, Ua] = theta_2_1_1(tinanew, tin, active, Hold, Ua, Uin, Wfdobj, Slambda, Hanew, N, tQin, crit);

end

theta_in(ind) = tin;
theta_out = theta_in;



