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

%% Begin theta_0_0_0()

[iter, tina, active, Ua, Hold] = theta_0_0_0(tin, Wfdobj, Uin, N, Slambda, tQin, crit);

% Wmat   = eval_fd(tin, Wfdobj);
% DWmat  = eval_fd(tin, Wfdobj, 1);
% D2Wmat = eval_fd(tin, Wfdobj, 2);
% Emat   = exp(Wmat);
% Pmat   = Emat./(1+Emat);
% Rmat   = Uin - Pmat;
% 
% %  initial function, first and second derivative values  wrt theta
% 
% Hold   = NaN*ones(Nin,1);
% DHold  = zeros(Nin,1);
% D2Hold = zeros(Nin,1);
% 
% Hold   = -sum(Uin.*Wmat - log(1+Emat),2)/N;
% DHold  = -sum(Rmat .*DWmat, 2)/N;
% D2Hold =  sum(-D2Wmat.*Rmat + DWmat.^2.*Pmat.*(1-Pmat), 2)/N;
% if Slambda > 0
%     Hold   = Hold   +    Slambda.*(tin - tQin).^2/N;
%     DHold  = DHold  + 2.*Slambda.*(tin - tQin)/N;
%     D2Hold = D2Hold + 2.*Slambda/N;
% end
% 
% %  initialize iterations
% 
% iter = 0;
% 
% %  find cases where gradient size exceeds criterion and further 
% %  optimization is required
% 
% active  = find(abs(DHold) > crit & D2Hold > 0);
% %nactive = length(active);
% %  select active theta's and data
% tina    = tin(active);
% Ua      = Uin(active,:);  
% %  update scores in active set by a Newton-Raphson step
% 
% %troubleind = [982,1649];

% End theta_0_0_0()

while ~isempty(active) && iter < itermax
    
    %% Begin theta_0_1_0()
    
    [iter, step] = theta_0_1_0(iter, tina, Wfdobj, N, Ua, Slambda, tQin, active);
    % disp(step)
    
%     iter = iter + 1;
%     %  compute probabilities and other quantities for active cases
%     Wmata   = eval_fd(tina, Wfdobj);
%     DWmata  = eval_fd(tina, Wfdobj, 1);
%     D2Wmata = eval_fd(tina, Wfdobj, 2);
%     Emata   = exp(Wmata);
%     Pmata   = Emata./(1+Emata);
%     Rmata   = Ua - Pmata;
%     DHa     = -sum(Rmata.*DWmata, 2)/N;
%     D2Ha    =  sum(-D2Wmata.*Rmata + DWmata.^2.*Pmata.*(1-Pmata), 2)/N;
%     %  add regularization terms if required
%     
%     tQactive = tQin(active);
%     
%     if Slambda > 0
%         DHa  = DHa  + 2.*Slambda.*(tina - tQin(active))/N;
%         D2Ha = D2Ha + 2.*Slambda/N;
%     end
%     
%     %  disp([DHa, D2Ha])
%     %  find an step size that reduces all function values in
%     %  the active set
%     %  compute negative of initial step for each theta
%     step = DHa./D2Ha;
%     step(step < -2) = -2; 
%     step(step >  2) =  2; 
%     
    % End theta_0_1_0()
    
%     ind39  = find(active == 39);
%     ind309 = find(active == 309);
%     if ~isempty(ind39)
%           disp(['Iteration ',num2str(iter),' theta(39) = ', ...
%           num2str([tin(39),DHa(ind39),D2Ha(ind39),step(ind39)])])
%     end
%     if ~isempty(ind309)
%           disp(['Iteration ',num2str(iter),' theta(309) = ', ...
%           num2str([tin(309),DHa(ind309),D2Ha(ind309),step(ind309)])])
%     end
    %  ensure that initial step does not go to boundaries
    
    %% Begin theta_0_1_1()
    
    [tinanew, Hanew, fac, stepiter, failindex] = theta_0_1_1(tina, step, scrmin, scrmax, tin, active, Hold);
    
%     stepind0 = find(tina - step <= scrmin);
%     stepindn = find(tina - step >= scrmax);
%     if ~isempty(stepind0)  
%         step(stepind0) =  (tina(stepind0) - scrmin*1.01);  
%     end
%     if ~isempty(stepindn) 
%         step(stepindn) = -(scrmax*0.99 - tina(stepindn));
%     end
%     %  half fac as required to achieve the reduction of all fn. values
%     %  step through the iteration while any(Hold(active) - Ha < -1e-2)
%     %  up to a maximum number of 10 step halvings
%     %  the slack -1e-2 was used to reduce the number of halvings, and
%     %  was chosen after experimentation with different values
%     tinanew  = tin(active);
%     Hanew    = Hold(active);
%     fac      = 1;
%     stepiter = 0;
%     %  initial set of theta's failing to meet criterion
%     failindex  = ones(nactive,1);
%     failindex  = logical(failindex);
    
    % End theta_0_1_1()
    
    while sum(failindex) > 0 && stepiter < 10
        
        %% Begin theta_0_2_0()
        
        [stepiter, tinanew] = theta_0_2_0(stepiter, step, scrmin, scrmax, failindex, fac, tinanew, active);
        
%         stepiter = stepiter + 1;
%         % disp(num2str([stepiter, fac]))
%         %  compute current step size
%         activefail  = active(failindex);
%         stepnew     = fac*step(failindex);
%         tinanewfail = tinanew(failindex);
%         stepnewind0 = find(tinanewfail - stepnew <= scrmin);
%         stepnewindn = find(tinanewfail - stepnew >= scrmax);
%         if ~isempty(stepnewind0)  
%             stepnew(stepnewind0) = ...
%                  (tinanewfail(stepnewind0) - scrmin*1.01);  
%         end
%         if ~isempty(stepnewindn)  
%             stepnew(stepnewindn) = ...
%                 -(scrmax*0.99 - tinanewfail(stepnewindn));  
%         end
%         %  update ability values
%         tinanew(failindex) = tinanew(failindex) - stepnew;
        
        % End theta_0_2_0()
        
        %% Begin theta_0_2_1()
        
        [Hanew, failindex, fac] = theta_0_2_1(Hanew, failindex, fac, Wfdobj, tinanew, Ua, Slambda, tQin, N, Hold, crit, active);
        
%         %  Compute function value
%         Wmata(activefail,:) = eval_fd(tinanew(failindex), Wfdobj);
%         wfail = Wmata(activefail,:);
%         Emata(activefail,:) = exp(Wmata(activefail,:));
%         %  current vector of function values
%         Hanew(failindex) = -sum(Ua(failindex,:).*Wmata(activefail,:) - ...
%                                 log(1+Emata(activefail,:)),2)/N;
%         if Slambda > 0
%             Hanew(failindex) = Hanew(failindex) + ...
%                 Slambda.*(tinanew(failindex) - ...
%                           tQin(activefail)).^2/N;
%         end
%         % disp(num2str([activefail, tinanew(failindex), ...
%         %               Hold(activefail),Hanew(failindex)]))
%         failindex = Hold(active) - Hanew < -crit;
%         %  halve fac in preparation for next step size iteration
%         fac = fac/2;
        
        % End theta_0_2_1()
    end
    
    %% Begin theta_0_3_0()
%     if stepiter == 10
%         warning('Number of halvings = 10.')
%     end

    [tina, tin, active, Hold, Ua] = theta_0_3_0(tinanew, tin, active, Hold, Ua, Uin, Wfdobj, Slambda, Hanew, N, tQin, crit);

    %  either all function values reduced or 10 halvings completed
    %  save current function values
%     tin(active) = tinanew;
%     Wmata   = eval_fd(tinanew, Wfdobj);
%     DWmata  = eval_fd(tinanew, Wfdobj, 1);
%     D2Wmata = eval_fd(tinanew, Wfdobj, 2);
%     Emata   = exp(Wmata);
%     Pmata   = Emata./(1+Emata);
%     Rmata   = Ua - Pmata;
%     DHanew  = -sum(Rmata.*DWmata, 2)/N;
%     D2Hanew =  sum(-D2Wmata.*Rmata + DWmata.^2.*Pmata.*(1-Pmata), 2)/N;
%     if Slambda > 0
%         DHanew  = DHanew  + 2.*Slambda.*(tinanew - tQin(active))/N;
%         D2Hanew = D2Hanew + 2.*Slambda/N;
%     end
%     Hold(active)   = Hanew;
%     DHold(active)  = DHanew;
%     D2Hold(active) = D2Hanew;
%     active  = find(abs(DHold) > crit & D2Hold > 0);
%     nactive = length(active);
%     if nactive > 0
%         tina = tin(active);
%         Ua   = Uin(active,:);
%     end
    %  return for a new optimization step
    
    % End theta_0_3_0()
end

% if iter == itermax
%     warning('Max. iterations reached, size of active follows:');
%     disp(active')
%     disp(tin(active)')
%     disp(Hold(active)')
%     disp(D2Hold(active)')
% %     hist(log10(D2Hold))
% %     pause
% end

%  opimization complete for all cases, save results and exit

theta_in(ind) = tin;
theta_out = theta_in;



