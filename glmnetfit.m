

function [Theta, I, LL, badcond,lgm, max_iterations_reached,design_is_singular, N] = mnlfit(R,Ysp,lambda,alpha,varargin)

%function [Theta, I] = mnlfit(R,Ysp,lambda,alpha'param',value,... )
%
%     Fits a multinomial generalized linear model with regressor R and
%     observation vector Ysp, using regularized glmnet procedure (See eg
%     Friedman, Hastie and Tibshirani 2010). Parameters are cyclically, one
%     at a time, cycling optionallyover parameters and optionally over classes.
%
%     R is a regressor structure as returned by makeregresso, where 
%     R.value is an N x M matrix.
%
%     Ysp is the Mx1 vector of observations.with 1 denoting that the chosen option
%     associated with the corresponding row of X.
%
%     lambda is strength of regularization and alpha is the trade off
%     between l2 (alpha=1) and l1 (alph=0) regularization.
%
%     In the case of multiple assignment, the outcome is
%     treated as a weighted mixture of the assigned values. In that case 
%     the value of Y(i) gives the relative weight accorded the corresponding option. 
%
%     SEE ALSO makeregressor, modelFit
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2010


InitTheta = [];
Hreg= 0;
Lreg = 0;
runiter=true;
fix = [];
firth = false;
debugstop = false;
checkdesign = false;
diagsonly = [];
surrogate= true;
LC = [];
binvolume = 0;
max_iterations_reached = 0;
design_is_singular = [];
LevMar = true; %use Levenberg-Marquardt damping if true
showprog = false; %Generate output % Fitting is ~6 x faster without display!
maxcount = 500; %maximum iterations to continue 
obsfreq = 1;
discard = [];
L2reg = 0;
showstep = 20;
cycle_classes = true;
% LC = 0;
i=1;

L1reg = 0;

L2reg = 0;

tol = 1e-6;

%%%% Options are described below %%%%
while i <= length(varargin)
    
    switch lower(varargin{i})
        case {'gaussreg','regularization'}   %Specify gaussian regularization  the SQUARE of the weighting on the L2 norm (sqrt(g)*|r|^2)
            Hreg = varargin{i+1};
            i = i+1;            
 
        case 'runiter'    %do not run iterations (ie only compute likelihood at InitTheta) if false
            runiter = varargin{i+1};
            i = i+1;
        case 'maxiter'    %maximum iterations
            maxcount = varargin{i+1};
            i = i+1;
        case 'inittheta'     %Only fits the full model (no LLR statistic on submodels
            InitTheta = varargin{i+1};
            i = i+1;
        case 'binvolume'     %Adds the specified vector of constants to log weight to correct for bin volume.
            binvolume = varargin{i+1};
            i = i+1;
        case 'fix'     %Leave specified parameters fixed (not subject to maximization); Change this to use constraint.
            fix = varargin{i+1};
            i = i+1;
        case 'obsfreq'     % Observation frequency: A vector with the frequency of each observation (the default is 1). 
                           % Collapsing over trials with equivalent design and outcome can greatly increase efficiency.  
            obsfreq = varargin{i+1};
            i = i+1;
       
        case 'diagsonly'%ignore cross terms in computing 2nd derivative
            diagsonly= varargin{i+1};
            i = i+1;         
        case 'surrogate'     %
            surrogate = varargin{i+1};
            i = i+1;
        case 'checkdesign'     % checks if the design matrix is singular
            checkdesign = varargin{i+1};
            i = i+1;
        case 'linearconstraint'     %a matrix of linear constraints on maximization, which is a contrast matrix such that
              LC = varargin{i+1};   %  a*C = 0 if c is the same length as
              i = i+1;              %  the parameters and a*C - C(end,:) = 0
                                    %  if it is longer by 1. Lagrange
                                    %  multipliers are appended to the
                                    % parameter vector.
        case 'include_null'     %Doesn't do anything now
            const = varargin{i+1}; %#ok<NASGU>
            i = i+1;
        case {'show_progress','showprog'} 
            showprog = varargin{i+1};
            i = i+1;
        case 'discard'     %
            discard = varargin{i+1};
            i = i+1;
        case 'showstep'     %
            showstep = varargin{i+1};
            i = i+1;
        case {'tolerance','tol'}     %
            tol = varargin{i+1};
            i = i+1;
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end
%%%%%%%%%%%%%%%%%%


L1reg = lambda*alpha;
L2reg = lambda*(1-alpha);

thresh=L1reg; % Soft threshold

if isempty(binvolume)
    binvolume = 0;
end

X = R.value;
b = R.noptions;
Npar = R.Npar;
% fixed = R.fixed;
if size(X,1) ~= length(Ysp)
    error('X must have the same number of rows as elements in Y')
end
    


% blocksum = sparseblock(ones(1,size(X,2)),b)';







if ~isempty(discard)
      blocksum = sparseblock(ones(1,size(X,1)),b)';
   disclarge = blocksum*discard;
   Ysp = Ysp(~disclarge);
   X = X(~disclarge,:);
   b = b(~discard); 
end

if length(obsfreq) > 1 
    if ~isempty(discard)
        obsfreq = obsfreq(~discard);
    end
    
      blocksum = sparseblock(ones(1,size(X,1)),b)';
    OF = diag(sparse(obsfreq)); %Frequency of each observation
    OFlarge = diag(sparse(blocksum*obsfreq)); %Frequency of each observation
    
%     L1reg = L1reg*sum(obsfreq);  %%% For lasso and ridge regression, prior is proportional to the sample size
%     L2reg = L2reg*sum(obsfreq);
else
    OF      = obsfreq(1);
    OFlarge = obsfreq(1);
    
%    L1reg = L1reg*length(b)*obsfreq(1);
%    L2reg = L2reg*length(b)*obsfreq(1);
    
end


Yblocksum = sparseblock(full(Ysp)',b)';
Y = full(Ysp);


clear R;

if firth && ( any(b>2) || any(any(abs(X(1:2:end,:).*X(2:2:end,:)>max((1e-10).*mean(X))) )) )
      error('Jeffrey''s prior is currently implemented only for dichotomous logistic regression.')
end

if checkdesign
    xx = X'*X;  
    design_is_singular = rank(xx) < length(xx);
    if design_is_singular
       beep
       fprintf('Design Matrix is singular or poorly conditioned! Press Ctrl-C to stop, or wait to continue')
       pause(2)
    end
end

%Create constraint matrix for fixed parameters
if any(fix)
%     constraint = true;
    fixed = find(fix);
    LC = zeros(Npar,1);
    for i = 1:length(fixed)
        LC(i,fixed(i)) = 1;
    end
        LC(end+1,:) = InitTheta(fixed);
end

transposeSpx = true;
% maxMatSize = 1e9./64; % Uses a slower loop to compute Hessian matrix if the number of non-zero elements exceeds this number

lmnu = 1.5; %Damping parameter 

badcond =0;

% % trblockstep = 1;

% 
% if nargin < 9 || isempty(surrogate) % Use a surrogate matrix in optimization rather than the actual 2nd derivative
%     surrogate = true;              % Fitting is faster
% end
recomputeD2Lintvl = 20; %If surrogate is true, 2nd derivative matrix is recomputed once every recomputeD2Lintvl iterations;
    

if nargin < 5 || isempty(runiter)
    runiter = true;
end

if  isempty(InitTheta) 
    InitTheta = zeros(Npar,1) + eps;
elseif length(InitTheta ) < Npar
    InitTheta(end+1:Npar) = 0;
end

if length(Hreg) == 1  %Gaussian prior
    Hreg = eye(Npar).*Hreg;
elseif isvector(Hreg)
    Hreg = diag(Hreg);
end
    

if length(L2reg) == 1  %Gaussian prior
    Hreg = Hreg + eye(Npar).*L2reg;
elseif isvector(L2reg)
    Hreg = Hreg + diag(L2reg);
else
    Hreg = Hreg + L2reg;
end

    

if length(Lreg) == 1  %Laplace prior
    Lreg = eye(Npar).*Lreg;
elseif isvector(Lreg)
    Lreg = diag(Lreg);
end



if length(L1reg) == 1  %L1 norm
    L1reg = ones(Npar,1).*L1reg;

end



if isempty(X)
    runiter = false;
    InitTheta = [];
end

if LevMar
    damp = .1;
    Slev = damp*eye(length(InitTheta));
else 
    damp = 0; %#ok<UNRCH>
end

if transposeSpx 
    X = X';
    binvolume = binvolume';
end    


Theta = InitTheta(:);

dstep = Inf;
count = 0;
if showprog
    fprintf('%5i: damping: %1.2e,  dstep = %1.2e',0,lmnu,0);
end

del = Inf;




%Second level for multiple assignment
b2 = full(sum(Yblocksum));
b2 = b2(b2>1);
b2bl = blocknorm(full(Ysp)',b);

X2 = X(:,Ysp'.*b2bl>1);

if firth
      Cxx = repmat(X,Npar,1).*kron(X,ones(Npar,1));
%       Cxx2 = repmat(X2,Npar,1).*kron(X2,ones(Npar,1));
end

% if ~diagsonly
spX = sparseblock(X,b);
spX2 = sparseblock(X2,b2);
% end
if firth %Firth's correction requires that the information matrix be computed at every iteraion anyway
    surrogate = false;
end


if ~isempty(LC)    
    if size(LC,1) == Npar
        LC(end+1,:) = 0;
    end
    if size(LC,1) ~= Npar+1
        error('Contraint matrix must have as many rows as the number of parameters')
    end
    constraint = 1;
    
    lgm = zeros(size(LC,2),1);
else
    lgm = []; %lagrange multipliers    
    constraint = 0;
end

nlgm = size(LC,2);

blockn = blocknorm(ones(1,size(X,2)),b);

%%% weighted length of a vector
wleng = @(x, M ) sqrt(sum(x'*M*x)+eps);

parindex=1:Npar;
subi = @(X,indx) X(parindex==indx,parindex==indx);
%%% A matrix which will be used to compute 0th and 1st  derivatives of
%%% the prior
RegMat = @(theta,a,indx) ( a*subi(Hreg,indx) + subi(Lreg,indx)./wleng(theta(indx),subi(Lreg,indx)) )*theta(indx); 

%%% Functions to compute log likelihood and outcome probability
lrrfun = @(theta) theta' * X + binvolume; %compute the log weight
nlrrfun = @(theta) lrrfun(theta) - blocknorm(lrrfun(theta),b)./blockn; %subtract mean to improve numerical stability
Pfun =  @(theta) exp( nlrrfun(theta) )./  (blocknorm(exp( nlrrfun(theta) ) ,b)); %Proabability of each outcome

LLfun = @(theta) sum(log(Pfun(theta)*Yblocksum)*OF) - theta'*RegMat(theta,1,parindex) - L1reg(:)'*abs(theta); %Log Likelihood

%%% Function to compute the second derivative of the prior distribution

if Lreg == 0
    D2Regfun =@(theta,indx) subi(Hreg,indx); %%% Gaussian prior only
else
    D2Regfun =@(theta,indx)  Lreg/wleng(theta(indx),subi(Lreg,indx)) - (subi(Lreg,indx)*theta(indx)*theta(indx)'*subi(Lreg,indx)')/wleng(theta(indx),subi(Lreg,indx))^3   + subi(Hreg,indx); %%% Gaussian and laplace prior 
end


if all(L1reg == 0) % apply l1 norm
    DLfun = @(theta,P,P2,indx)  X(indx,:)*OFlarge*(P2 - P)' - RegMat(theta,2,indx);
else
    DLfun = @(theta,P,P2,indx)  X(indx,:)*OFlarge*(P2 - P)' - RegMat(theta,2,indx) - L1reg(indx).*sign(theta(indx));
end    

stop = false;
lastlap = false;


while  ~stop  && runiter  
    
    Theta0=Theta;
    for k = 1:Npar
        P =  Pfun(Theta); %Proabability of each outcome

        W = P.*(1-P);

        Pnorm = blocknorm(P.*Y',b);
        P2 = P.*Y'./Pnorm;


        %%% Derivative of the log likelihood

    %     DL = X*OFlarge*(P2 - P)' - RegMat(Theta,2);
        DL = DLfun(Theta,P,P2,k);


        if constraint 
            DL = DL + LC(1:end-1,:)*lgm; %Derivative of parameters under linear constraint
            Dlgm = [Theta; -1]'*LC;             % Derivative of Lagrange Multipliers;
        end


            if debugstop && (isnan(LL1) || isnan(LL2))
                keyboard
            end
            
        S2 = 0;
        X2 = X(:,b2bl.*Y' > 1);
%          P2 = P2(:,b2bl.*Y' > 1);



        if ~isempty(X2) 

                dgP2 = sparseblockmex(P2,0:length(P2), 0:length(P2),length(P2)); %Equivalent to but faster than diag(P) or diag(sparse(P))
                spP2 = sparseblock(P2,b2);
                S2 = full( X2*dgP2*X2'  - spblocktrace( (spX2*spP2')*(spP2*spX2'), Npar ));
       
        end

           dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
           spP = sparseblock(P,b);
           spX = sparseblock(X(k,:),b);
           
%            S1 = full( - X(k,:)*dgP*OFlarge*X(k,:)'  + spblocktrace( (spX*spP')*OF*(spP*spX'), Npar));
           S1 = full( - X(k,:)*dgP*OFlarge*X(k,:)'  + trace( (spX*spP')*OF*(spP*spX')));
            D2L = S2 + S1 - D2Regfun(Theta,k);        
           dth=-X(k,:)*(P2-P)'./D2L;
           Theta(k)=Theta(k)+dth;
           
           %%% Soft thresholding..
           Theta(k)=Theta(k)*(abs(Theta(k))>L1reg(k));
           
%        dstep = abs(dll);

    end
      de=sqrt( damp*sum(del.^2)./sum(Theta.^2));
  if showprog && mod(count,showstep) == 0
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
        if ispc 
            fprintf('\b\b') %%% For some reason, needs extra backspaces on PCs            
        end
        if ~isinf(dstep)
            fprintf('%5i: damping: %1.2e,  dstep = %1.2e',count,0,de);
        else
            fprintf('%5i: damping: %1.2e,  dstep = %1.2e',count,0,0);
        end
    end        

    count = count+1;
    
%     dll = LLfun(Theta)-LLfun(Theta0);

%     de=sqrt( damp*sum(del.^2)./sum(Theta.^2));
     if count > maxcount % && std(llhist(end-50:end))./diff(llhist([1,end])) < 1e6 

         warning('GazeReader:NoConvergence','Failed to converge after %i iterations. Estimate may be unbounded.',maxcount)
         pause(1)
         max_iterations_reached = 1;
         break
     end

    del = Theta-Theta0;
    dstep = 0;
    
     %%% avoid stopping too early if the Hessian needs to be updated. 
%      if dstep + sqrt( damp*sum(del.^2)./sum(Theta.^2)) < tol %%% damping dependent term avoids stopping when the damping becomes high
     if dstep + sqrt( damp*sum(del.^2)./sum(Theta.^2)) < tol %%% damping dependent term avoids stopping when the damping becomes high
         if lastlap
             stop = true;
         else
             lastlap = true;
         end         
     else
         lastlap = false;
     end

end

if showprog, fprintf('\n');end

%One more pass to Compute final value of everyting and observed Fisher information 
% lrr = Theta' * X + binvolume;
% lrr = lrr - blocknorm(lrr,b)./blockn; %subtract mean to improve numerical stability
% P =  exp( lrr )./  (blocknorm(exp( lrr ) ,b)); %Proabability of each outcome
% LL = sum(log(P*Yblocksum+eps)*OF);
% LL = sum(log(P*Yblocksum+eps)*OF) - Theta'*RegMat*Theta;

P = Pfun(Theta);
LL = LLfun(Theta);

Pnorm = blocknorm(P.*Y',b);
P2 = P.*Y'./Pnorm;
% DL = X*OFlarge*(P2 - P)' - RegMat(Theta,2);

X2 = X(:,b2bl.*Y' > 1);
P2 = P2(:,b2bl.*Y' > 1);
dgP2 = sparseblockmex(P2,0:length(P2), 0:length(P2),length(P2)); %Equivalent to but faster than diag(P) or diag(sparse(P))
spP2 = sparseblock(P2,b2);
S2 = full( X2*dgP2*X2'  - spblocktrace( (spX2*spP2')*(spP2*spX2'), Npar  ));
dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
spP = sparseblock(P,b);
D2L = full( S2 - X*dgP*OFlarge*X'  + spblocktrace( (spX*spP')*OF*(spP*spX'), Npar ))- D2Regfun(Theta);


if firth
    LL = LL + .5*log(det(-D2L));         
end

N = sum(ones(1,length(b))*OF);

if any(P==1) || any(P==0)
            fprintf('\nProbability is saturated. Model may be unbounded.\n')
end
       
if runiter
    I = -D2L;             
    if any(isnan(I(:))) || rank(I) < length(I)
        badcond = 1;
    end

else
    I = nan;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%

