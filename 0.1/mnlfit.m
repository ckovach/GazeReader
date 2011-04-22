

function [Theta, I, LL, badcond,lgm, max_iterations_reached,design_is_singular, N] = mnlfit(R,Ysp,varargin)

%function [Theta, I] = mnlfit(R,Ysp,'param',value,... )
%
%     Fits a multinomial generalized linear model with regressor R and
%     observation vector Ysp, using damped Newton-Raphson gradient ascent
%     on the likelihood function.
%
%     R is a regressor structure as returned by makeregresso, where 
%     R.value is an N x M matrix.
%
%     Ysp is the Mx1 vector of observations.with 1 denoting that the chosen option
%     associated with the corresponding row of X.
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

% C. Kovach 2008


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
showprog = true; %Generate output 
maxcount = 500; %maximum iterations to continue 
obsfreq = 1;
discard = [];
% LC = 0;
i=1;

%%%% Options are described below %%%%
while i <= length(varargin)
    
    switch lower(varargin{i})
        case {'gaussreg','regularization'}   %Specify gaussian regularization (ridge)
            Hreg = varargin{i+1};
            i = i+1;            
        case {'laplreg'}   %Specify laplacian regularization
            Lreg = varargin{i+1};
            i = i+1;            
        case 'firth'    %Use Jeffrey's prior as described by Firth: currently 
                        % works only for logisitic regresion normalized to 2nd option!
            firth = varargin{i+1};
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
        case 'show_progress' 
            showprog = varargin{i+1};
            i = i+1;
        case 'discard'     %
            discard = varargin{i+1};
            i = i+1;
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end
%%%%%%%%%%%%%%%%%%

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
else
    OF      = obsfreq(1);
    OFlarge = obsfreq(1);
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
       pause(5)
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
    

if length(Lreg) == 1  %Laplace prior
    Lreg = eye(Npar).*Lreg;
elseif isvector(Lreg)
    Lreg = diag(Lreg);
end

% if any(fixed~=0)      %Check Regressor structure to see if any values are to remain fixed.
%     fix = fixed~=0;
%     InitTheta(fix) = fixed(fix);
% end

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

tol = 1e-6;

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

RegMat = @(theta,a) ( a*Hreg + Lreg./sqrt(sum(theta.^2)+eps) )*theta;  %Compute a matrix to adjust 1st and 2nd derivatives by the prior distirbution.

%%% Functions to compute log likelihood and outcome probability
lrrfun = @(theta) theta' * X + binvolume; %compute the log weight
nlrrfun = @(theta) lrrfun(theta) - blocknorm(lrrfun(theta),b)./blockn; %subtract mean to improve numerical stability
Pfun =  @(theta) exp( nlrrfun(theta) )./  (blocknorm(exp( nlrrfun(theta) ) ,b)); %Proabability of each outcome
LLfun = @(theta) sum(log(Pfun(theta)*Yblocksum)*OF) - theta'*RegMat(theta,1); %Log Likelihood

while dstep + sqrt( damp*sum(del.^2)./sum(Theta.^2)) > tol  && runiter   ; %Added damping dependent term to avoid premature halting when damping is high
% while dstep > tol  && runiter   

%     lrr = Theta' * X + binvolume;
%     lrr = lrr - blocknorm(lrr,b)./blockn; %subtract mean to improve numerical stability
%     P =  exp( lrr )./  (blocknorm(exp( lrr ) ,b)); %Proabability of each outcome

    P =  Pfun(Theta); %Proabability of each outcome
    
%     P =  exp( Theta' * X + binvolume)./  (blocknorm(exp( Theta' * X + binvolume ) ,b)+eps); %Proabability of each outcome

    Pnorm = blocknorm(P.*Y',b);
    P2 = P.*Y'./Pnorm;
    
    % Derivative of the log likelihood
    
    DL = X*OFlarge*(P2 - P)' - RegMat(Theta,2);
    
    if constraint 
        DL = DL + LC(1:end-1,:)*lgm; %Derivative of parameters under linear constraint
        Dlgm = [Theta; -1]'*LC;             % Derivative of Lagrange Multipliers;
    end
    
    if ~surrogate || count == 0 || mod(count,recomputeD2Lintvl)==0
        
         X2 = X(:,b2bl.*Y' > 1);
         P2 = P2(:,b2bl.*Y' > 1);


        S2 = 0;
        if ~isempty(P2) 

            if diagsonly
                S2 = X2*diag(sparse(P2.*(1-P2)))*X2';
            else
                dgP2 = sparseblockmex(P2,0:length(P2), 0:length(P2),length(P2)); %Equivalent to but faster than diag(P) or diag(sparse(P))
                spP2 = sparseblock(P2,b2);
                S2 = full( X2*dgP2*X2'  - spblocktrace( (spX2*spP2')*(spP2*spX2'), Npar ));
            end

        end

        % Computing the second derivative

        if diagsonly 
            
            W = sparseblockmex(P.*(1-P),0:length(P), 0:length(P),length(P))*OFlarge; %Equivalent to but faster than diag(P) or diag(sparse(P))
            S1 = - X*W*X';
            D2L = S2 + S1 - Hreg;        
        else 
           dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
           spP = sparseblock(P,b);
           
%            q = spX*spP'; %faster but uses more memory
%            S1 = full( - X*dgP*X'  + spblocktrace( q*q', Npar ));
%            clear q;
           S1 = full( - X*dgP*OFlarge*X'  + spblocktrace( (spX*spP')*OF*(spP*spX'), Npar));
           D2L = S2 + S1 - Hreg;        

        end
%         D2Linv = (D2L-Hreg)^-1;
    end
    
    if firth  %This has very limited functionality in the current version.

        %Modified to account for multiple assignment. This needs to be
        %double checked!! This still needs to be modified to accomodate
        %constraints!

%         Iinv = -(D2L - Hreg)^-1; 
        Iinv1 = -(S1 - Hreg)^-1; 
        
        
        DI1 = Cxx*(repmat((1-2*P).*P.*(1-P),Npar,1).*X)';
        
%         if ~isempty(P2)
%             Iinv2 = -(S2)^-1; 
%             DI2 = Cxx2*(repmat((1-2*P2).*P2.*(1-P2),Npar,1).*X2)';
%         else
%             Iinv2 = 0;
%             DI2 = 0;
%         end
        
        % For now, ignoring the second term in I related to multiple assignment
        DL = DL + .5*DI1'*Iinv1(:) ;% + .5*DI2'*Iinv2(:); %Modified score function by Jeffreys prior (see Firth 1993)
        
%         dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
%         spP = sparseblock(P,b);
%         S1 = full( - X*dgP*X'  + spblocktrace( (spX*spP')*(spP*spX'), Npar ));
%         D2L =  S1; % + S2;        

%         LL = sum(log(P*Yblocksum)*OF) + .5*log(det(-D2L));         
        LL = LLfun(Theta) + .5*log(det(-D2L));     
    else
        
        %Redundant. 
%         lrr = Theta' * X + binvolume;
%         lrr = lrr - blocknorm(lrr,b)./blockn; %subtract mean to improve numerical stability
%         P =  exp( lrr )./  (blocknorm(exp( lrr ) ,b)); %Proabability of each outcome

%         LL = sum(log(P*Yblocksum)*OF);
        LL = LLfun(Theta);
    end
    
    % Newton's method
    if ~LevMar  
        if ~constraint
%             del = -D2L^-1 * (DL - RegMat(Theta,2) ); %Hreg and Lreg are Gaussian and Laplacian priors, respectively.
            del = -D2L^-1 * DL; %Hreg and Lreg are Gaussian and Laplacian priors, respectively.
        else
            %Updating with lagrange multiplier
            Q = cat(2,D2L,LC(1:end-1,:));
            Q = cat(1,Q,[LC(1:end-1)',zeros(nlgm)]);
            
%             del = -Q^-1 * ( cat(1,DL,Dlgm) - cat(RegMat(Theta,2),zeros(nlgm,1)));
            del = -Q^-1 *  cat(1,DL,Dlgm) ;
            dellgm = del((1:nlgm)+Npar);
            del = del(1:Npar);
%             del(fix) = 0;
            
            lgm = lgm + dellgm;
        end            
        Theta = Theta + del;  %New Theta

    else    %if ~surrogate
        
        if ~constraint
%             del1 = -(D2L - Slev)^-1  * (DL - RegMat(Theta,2));
%             del2 = -(D2L - Slev./lmnu)^-1 * (DL - RegMat(Theta,2));
            del1 = -(D2L - Slev)^-1  * DL ;
            del2 = -(D2L - Slev./lmnu)^-1 * DL;
            
            Theta1 = Theta + del1;
            Theta2 = Theta + del2;
            
        else
            
            Qfun = @(nm) cat(1, cat( 2 , D2L - Slev/nm ,LC(1:end-1,:) ),...
                       [ LC(1:end-1,:)' , zeros(nlgm) ]  );
                   
%             Q1 = cat(2,D2L - Slev,LC(1:end-1,:));
%             Q1 = cat(1,Q1,[LC(1:end-1,:)',zeros(nlgm)]);
%             Q2 = cat(2,D2L - Slev./lmnu,LC(1:end-1,:));
%             Q2 = cat(1,Q2,[LC(1:end-1,:)',zeros(nlgm)]);

%             Q1 = Qfun(1);
%             Q2 = Qfun(lmnu);
            
%             delfun = @(mn) -Qfun(mn)^-1 * ( cat(1,DL,Dlgm(:)) - cat(1,RegMat(Theta,2),zeros(nlgm,1)));
            delfun = @(mn) -Qfun(mn)^-1 *  cat(1,DL,Dlgm(:));
            
%             del1 = -Q1^-1 * ( cat(1,DL,Dlgm(:)) - cat(1,RegMat(Theta,2),zeros(nlgm,1)));
            del1 = delfun(1);
            dellgm1 = del1((1:nlgm)+Npar);
            del1 = del1(1:Npar);

            
%             del2 = -Q2^-1 * ( cat(1,DL,Dlgm(:)) - cat(1,RegMat(Theta,2),zeros(nlgm,1)));
            del2 = delfun(lmnu);          
            dellgm2 = del2((1:nlgm)+Npar);
            del2 = del2(1:Npar);
            
            Theta1 = Theta + del1;
            Theta2 = Theta + del2;
            
        end            

        if ~firth
%             lrr = Theta1' * X + binvolume;
%             lrr = lrr - blocknorm(lrr,b)./blockn; %subtract mean to improve numerical stability
            
        
%             P =  exp( lrr )./  (blocknorm(exp( lrr ) ,b)); %Proabability of each outcome
%             LL1 = sum(log(P*Yblocksum)*OF) - Theta1'*RegMat(Theta1,1);

            LL1 = LLfun(Theta1);
            LL2 = LLfun(Theta2);

%            lrr = Theta2' * X + binvolume;
%             lrr = lrr - blocknorm(lrr,b)./blockn; %subtract mean to improve numerical stability
%             P =  exp( lrr )./  (blocknorm(exp( lrr ) ,b)); %Proabability of each outcome
% %      P =  exp( Theta2' * X + binvolume )./  (blocknorm(exp( Theta2' * X + binvolume) ,b) + eps); 
% %             LL2  = sum(log(P*Yblocksum)*OF);
%             LL2  = sum(log(P*Yblocksum)*OF) - Theta2'*RegMat(Theta2);

        else

            lrr = Theta1' * X + binvolume;
            lrr = lrr - blocknorm(lrr,b)./blockn; %subtract mean to improve numerical stability
            P =  exp( lrr )./  (blocknorm(exp( lrr ) ,b)); %Proabability of each outcome
%             P =  exp( Theta1' * X + binvolume)./  blocknorm(exp( Theta1' * X + binvolume) ,b); 
            dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
            spP = sparseblock(P,b);
            S1 = full( - X*dgP*OFlarge*X'  + spblocktrace( (spX*spP')*OF*(spP*spX'), Npar ));
            D2L =  S1; % + S2;        

            LL1 = sum(log(P*Yblocksum)*OF) + .5*log(det(-D2L)) - Theta1'*RegMat(Theta1,2);                     

            lrr = Theta2' * X + binvolume;
            lrr = lrr - blocknorm(lrr,b)./blockn; %subtract mean to improve numerical stability
            P =  exp( lrr )./  (blocknorm(exp( lrr ) ,b)); %Proabability of each outcome
%         P =  exp( Theta2' * X + binvolume)./  blocknorm(exp( Theta2' * X + binvolume) ,b); 
            dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
            spP = sparseblock(P,b);
            S1 = full( - X*dgP*OFlarge*X'  + spblocktrace( (spX*spP')*OF*(spP*spX'), Npar));
            D2L =  S1; % + S2;        

            
            LL2 = sum(log(P*Yblocksum)*OF) + .5*log(det(-D2L)) - Theta2'*RegMat(Theta2,2);                     

        end            

        if debugstop && (isnan(LL1) || isnan(LL2))
            keyboard
        end

        if (LL > LL1 || isnan(LL1)) && ( LL > LL2 || isnan(LL2))
             damp = damp*lmnu;
 
        elseif LL1 > LL2 || isnan(LL2)
            
            Theta = Theta1;
            del = del1;
            if constraint
                lgm = dellgm1+lgm;
            end
   
        else
            Theta = Theta2;
%             TH = TH2;
            del = del2;
            damp = damp./lmnu;
            if constraint
                lgm = dellgm2+lgm;
            end
    
  
        end
        Slev = damp*eye(length(Theta));
        
        
    end

        
    
    dstep = sqrt( sum(del.^2)./sum(Theta.^2));
    
    if showprog
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
        if ~isinf(dstep)
            fprintf('%5i: damping: %1.2e,  dstep = %1.2e',count,damp,dstep);
        else
            fprintf('%5i: damping: %1.2e,  dstep = %1.2e',count,damp,0);
        end
    end        
        count = count+1;

%          dsteps(count) = dstep;         
%          llhist(count) = LL;        
%          Thetas(:,count) = Theta;

%          if length(llhist)> 1 && LL-llhist(end-1) <0, keyboard; end
         
             
         if count > maxcount % && std(llhist(end-50:end))./diff(llhist([1,end])) < 1e6 
            
             warning('GazeReader:NoConvergence','Failed to converge after %i iterations. Estimate may be unbounded.',maxcount)
             pause(1)
             max_iterations_reached = 1;
%              badcond = 1;
             break
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
X2 = X(:,b2bl.*Y' > 1);
P2 = P2(:,b2bl.*Y' > 1);
dgP2 = sparseblockmex(P2,0:length(P2), 0:length(P2),length(P2)); %Equivalent to but faster than diag(P) or diag(sparse(P))
spP2 = sparseblock(P2,b2);
S2 = full( X2*dgP2*X2'  - spblocktrace( (spX2*spP2')*(spP2*spX2'), Npar  ));
dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
spP = sparseblock(P,b);
D2L = full( S2 - X*dgP*OFlarge*X'  + spblocktrace( (spX*spP')*OF*(spP*spX'), Npar ))- Hreg;


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

