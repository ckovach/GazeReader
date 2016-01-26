

function [Theta, I, LL, badcond] = logitmodelsp(R,Ysp,InitTheta, Sreg, runiter, fix, firth )

%function [Theta, I] = logitmodelsp(X,Ysp,InitTheta, Sreg )
%
%     X is a block diagonal matrix. Each block, Bi is an b(i) x M matrix. M
%     is the number of regressors, which is fixed across trials. b(i) is the
%     number of possible outcomes on trial i, which varies from trial to trial.
%     Ysp is a sparse vector conatining 1 for each row in X which corresponds 
%     to an observed outcome.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2008



X = R.value;
b = R.noptions;
Npar = R.Npar;

clear R;

LevMar = true; %use Levenberg-Marquardt damping if true

transposeSpx = true;
% maxMatSize = 1e9./64; % Uses a slower loop to compute Hessian matrix if the number of non-zero elements exceeds this number

lmnu = 1.5; %Damping parameter 

badcond =0;

maxcount = 150; %maximum iterations to continue with non-increasing likelihood
% % trblockstep = 1;


if nargin < 7 || isempty(firth) %Bias correction using Firth's method
    firth = false;
end

if nargin < 6
    fix = [];
end


if nargin < 5 || isempty(runiter)
    runiter = true;
end

if nargin < 3 || isempty(InitTheta) || all(InitTheta == 0)
    InitTheta = zeros(Npar,1) + eps;
elseif length(InitTheta ) < Npar
    InitTheta(end+1:Npar) = 0;
end

if nargin < 4 || isempty(Sreg)
    Sreg = 0; %Regularization term
elseif length(Sreg) == 1
    Sreg = eye(Npar).*Sreg;
end

if isempty(X)
    runiter = false;
    InitTheta = [];
end

if LevMar
    damp = .1;
    Slev = damp*eye(length(InitTheta));

end

if transposeSpx 
    X = X';
end    

tol = 1e-6;

Theta = InitTheta(:);

dstep = Inf;
count = 0;
fprintf('%5i: damping: %1.2e,  dstep = %1.2e     ',0,lmnu,Inf);


del = Inf;

YX = X(:,Ysp==1);

blocksum = sparseblock(ones(1,size(X,2)),b)';
Yblocksum = sparseblock(full(Ysp)',b)';

P =  exp( Theta' * X )./  blocknorm(exp( Theta' * X ) ,b); %Proabability of each outcome
Pobs = P*Yblocksum;    %Probability of observed outcome

if firth
      Cxx = repmat(X,Npar,1).*kron(X,ones(Npar,1));
      LL = -Inf;
else
    LL = sum(log(Pobs));
end
while dstep > tol && runiter

    P =  exp( Theta' * X )./  blocknorm(exp( Theta' * X ) ,b); %Proabability of each outcome

    Pobs = P*Yblocksum;    %Probability of observed outcome

    %Second level for multiple assignment
    b2 = sum(Yblocksum);
    Pnorm = blocknorm(Pobs,b2 );
    P2 = P(Ysp==1)./Pnorm;
    
    % Derivative of the log likelihood
    
    DL = YX*P2' - X*P';
    
    
    X2 = YX(:,P2~=1);
    P2 = P2(:,P2~=1);
    
    if ~isempty(P2)
        S2 = -X2*diag(sparse(P2.*(1-P2)))*X2';
    else 
        S2 = 0;
    end
    
    % Computing the second derivative
    

    W = sparseblockmex(P.*(1-P),0:length(P), 0:length(P),length(P)); %faster than diag(P) or diag(sparse(P))
%     W2 = sparseblockmex(sqrt(P.*(1-P)),0:length(P), 0:length(P),length(P)); %faster than diag(P) or diag(sparse(P))

  
    D2L = S2 - X*W*X';
    
    if firth
        %This does not account for multople categorization yet!!!
%         H = dgP.*sqrt(1-dgP.^2)*X*(X*dgP.^2.*(1-dgP.^2)*X')^-1;
        Iinv = -D2L^-1; 
        
        
        DI = Cxx*(repmat((1-2*P).*P.*(1-P),Npar,1).*X)';
                
        DL = DL + .5*DI'*Iinv(:); %Modified score function by Jeffreys prior (see Firth 1993)
        LL = sum(log(P*Yblocksum)) + .5*log(det(-D2L));         
    else
        
        P =  exp( Theta' * X )./  blocknorm(exp( Theta' * X ) ,b); 
        LL = sum(log(P*Yblocksum));
    end
    
    % Newton's method
    
    if ~LevMar  
        del = -(D2L - Sreg)^-1 * (DL - Sreg*Theta);
        del(fix) = 0;

        Theta = Theta + del;
            % Computing the log likelihood
%         TH = repmat(diag(Theta),1,nblocks);

    
    else
        
%         del1 = -(D2L - Sreg*lmnu)^-1 * DL;
        del1 = -(D2L - Slev - Sreg)^-1  * (DL - Sreg*Theta);
        del1(fix) = 0;
        del2 = -(D2L - Slev./lmnu - Sreg)^-1 * (DL - Sreg*Theta);
        del2(fix) = 0;
        Theta1 = Theta + del1;
        Theta2 = Theta + del2;
        
        if ~firth
            P =  exp( Theta1' * X )./  blocknorm(exp( Theta1' * X ) ,b); 
            LL1 = sum(log(P*Yblocksum));
            P =  exp( Theta2' * X )./  blocknorm(exp( Theta2' * X ) ,b); 
            LL2  = sum(log(P*Yblocksum));
        else
            
            P =  exp( Theta1' * X )./  blocknorm(exp( Theta1' * X ) ,b); 
            I1 = -D2L; I2 = -D2L; 
            I1(:) = I1(:) + DI*del1; %First order approximation to new I's
            I2(:) = I2(:) + DI*del2; %This might affect convergence.
            LL1 = sum(log(P*Yblocksum)) + .5*log(det(I1)); 
            P =  exp( Theta2' * X )./  blocknorm(exp( Theta2' * X ) ,b); 
            LL2  = sum(log(P*Yblocksum)) + .5*log(det(I2));
        end            


        if LL > LL1 && LL > LL2
             damp = damp*lmnu;
 
        elseif LL1 > LL2
            
            Theta = Theta1;
            del = del1;
   
        else
            Theta = Theta2;
%             TH = TH2;
            del = del2;
            damp = damp./lmnu;
    
  
        end
        Slev = damp*eye(length(Theta));
    end
    dstep = sqrt( sum(del.^2)./sum(Theta.^2));
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    if ~isinf(dstep)
        fprintf('%5i: damping: %1.2e,  dstep = %1.2e',count,damp,dstep);
    else
        fprintf('%5i: damping: %1.2e,  dstep = %1.2e     ',count,damp,dstep);
    end
        
        count = count+1;

         dsteps(count) = dstep;
         
         llhist(count) = LL;
         
         if count > maxcount && std(llhist(end-50:end))./diff(llhist([1,end])) < 1e6 
            
             warning('Failed to converge after %i ieration. Estimate is likely unbounded.',maxcount)
             pause(1)
             badcond = 1;
             break
         end
end
fprintf('\n');


if runiter
    I = -D2L;            
else
    I = nan;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = makeSEtrial(b,Npar)

%Makes a diagonal block sparse matrix matrix of ones with each block b(i) x b(i)

Ncol = sum(b)';
mindex = zeros(Ncol*Npar,1);
mindex(Npar*cumsum(b(1:end-1))) = b(1:end-1);

modvec = zeros(Ncol*Npar,1);
modvec(1) = b(1);
modvec(mindex~=0) = diff(b);
modvec = cumsum(modvec);

strow = cumsum(mindex);

diter = ones(Npar*sum(b),1);
diter(1) = 0;
diter(Npar*cumsum(b(1:end-1))) =  1 - Npar*b(1:end-1);

iter = cumsum(diter);

Ir = mod(iter,modvec) + strow; 

Jc = find(mod(iter,modvec) == 0)-1;
Jc(end+1) = Npar*sum(b);

nrows = sum(b);

S = sparseblockmex(ones(sum(b.^2),1),Ir,Jc,nrows);


