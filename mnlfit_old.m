

function [Theta, I, LL, badcond] = mnlfit(R,Ysp,InitTheta, Sreg, runiter, fix, firth ,diagsonly,surrogate)

%function [Theta, I] = mnlfit(R,Ysp,InitTheta, Sreg )
%
%     R is a regressor structure as returned by makeregressor. 
%     R.value is an N x M matrix.
%     Ysp is the Mx1 vector of observations.with 1 denoting that the option
%     associated with the corresponding row of X was chosen.
%
%     In the case of multiple assignment, the outcome is
%     treated as a weighted combination of the assigned values. In that case 
%     the value of Y(i) gives the relative weight accoded the corresponding option. 

% C. Kovach 2008



X = R.value;
b = R.noptions;
Npar = R.Npar;
fixed = R.fixed;
if size(X,1) ~= length(Ysp)
    error('X must have the same number of rows as elements in Y')
end
    
clear R;

LevMar = true; %use Levenberg-Marquardt damping if true

transposeSpx = true;
% maxMatSize = 1e9./64; % Uses a slower loop to compute Hessian matrix if the number of non-zero elements exceeds this number

lmnu = 1.5; %Damping parameter 

badcond =0;

maxcount = 300; %maximum iterations to continue 
% % trblockstep = 1;


if nargin < 9 || isempty(surrogate) % Use a surrogate matrix in optimization rather than the actual 2nd derivative
    surrogate = true;              % Fitting is faster
end
recomputeD2Lintvl = 20; %If surrogate is true, 2nd derivative matrix is recomputed once every recomputeD2Lintvl iterations;

if nargin < 8 || isempty(diagsonly) % Ignore cross terms. 
    diagsonly = false;              % Converegence is slower but each iteration faster
%     diagsonly = true;              
end

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

if any(fixed~=0)      %Check Regressor structure to see if any values are to remain fixed.
    fix = fixed~=0;
    InitTheta(fix) = fixed(fix);
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
fprintf('%5i: damping: %1.2e,  dstep = %1.2e',0,lmnu,0);


del = Inf;


% blocksum = sparseblock(ones(1,size(X,2)),b)';
Yblocksum = sparseblock(full(Ysp)',b)';

% P =  exp( Theta' * X )./  blocknorm(exp( Theta' * X ) ,b); %Proabability of each outcome
% Pobs = P*Yblocksum;    %Probability of observed outcome


Y = full(Ysp);

%Second level for multiple assignment
b2 = full(sum(Yblocksum));
b2 = b2(b2>1);
b2bl = blocknorm(full(Ysp)',b);

X2 = X(:,Ysp'.*b2bl>1);

if firth
      Cxx = repmat(X,Npar,1).*kron(X,ones(Npar,1));
      Cxx2 = repmat(X2,Npar,1).*kron(X2,ones(Npar,1));
end

% if ~diagsonly
spX = sparseblock(X,b);
spX2 = sparseblock(X2,b2);
% end
if firth %Firth's correction requires that the information matrix be computed at every iteraion anyway
    surrogate = false;
end

% if surrogate % This computes a surrogate 2nd derivative matrix which is a lower bound of the actual one. 
%     P =  exp( zeros(size(Theta))' * X )./  blocknorm(exp( zeros(size(Theta))' * X ) ,b); 
%     Pnorm = blocknorm(P.*Y',b);
%     P2 = P.*Y'./Pnorm;
%     X2 = X(:,b2bl.*Y' > 1);
%     P2 = P2(:,b2bl.*Y' > 1);
%     dgP2 = sparseblockmex(P2,0:length(P2), 0:length(P2),length(P2)); %Equivalent to but faster than diag(P) or diag(sparse(P))
%     spP2 = sparseblock(P2,b2);
%     S2 = full( X2*dgP2*X2'  - spblocktrace( (spX2*spP2')*(spP2*spX2'), Npar ));
%     dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
%     spP = sparseblock(P,b);
%     D2L = (full( S2 - X*dgP*X'  + spblocktrace( (spX*spP')*(spP*spX'), Npar ) + Sreg));
% %     D2Linv = Hsurr^-1;
% end
% fprintf('%5i: damping: %1.2e,  dstep = %1.2e',0,0,0);

% while (dstep > tol || damp > 1e-1) && runiter 
while dstep > tol  && runiter 

    P =  exp( Theta' * X )./  blocknorm(exp( Theta' * X ) ,b); %Proabability of each outcome

    Pnorm = blocknorm(P.*Y',b);
    P2 = P.*Y'./Pnorm;
    
    % Derivative of the log likelihood
    
    DL = X*(P2 - P)';

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
            
            W = sparseblockmex(P.*(1-P),0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))
            S1 = - X*W*X';
            D2L = S2 + S1;        
        else 
           dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
           spP = sparseblock(P,b);

           S1 = full( - X*dgP*X'  + spblocktrace( (spX*spP')*(spP*spX'), Npar ));
           D2L = S2 + S1;        

        end
%         D2Linv = (D2L-Sreg)^-1;
    end
    
    if firth

        %Modified to account for multiple assignment. This needs to be
        %double checked!!

%         Iinv = -(D2L - Sreg)^-1; 
        Iinv1 = -(S1 - Sreg)^-1; 
        
        
        DI1 = Cxx*(repmat((1-2*P).*P.*(1-P),Npar,1).*X)';
        
%         if ~isempty(P2)
%             Iinv2 = -(S2)^-1; 
%             DI2 = Cxx2*(repmat((1-2*P2).*P2.*(1-P2),Npar,1).*X2)';
%         else
%             Iinv2 = 0;
%             DI2 = 0;
%         end
        
        % For now, ignoring the second term in I related to multiple assignment
        DL = DL - .5*DI1'*Iinv1(:) ;% + .5*DI2'*Iinv2(:); %Modified score function by Jeffreys prior (see Firth 1993)
        
        LL = sum(log(P*Yblocksum)) + .5*log(det(-D2L));         
    else
        
        P =  exp( Theta' * X )./  blocknorm(exp( Theta' * X ) ,b); 
        LL = sum(log(P*Yblocksum));
    end
    
    % Newton's method
    
    if ~LevMar  
        del = -D2L^-1 * (DL - Sreg*Theta);
        del(fix) = 0;

        Theta = Theta + del;

    else    %if ~surrogate
        
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
            dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
            spP = sparseblock(P,b);
            S1 = full( - X*dgP*X'  + spblocktrace( (spX*spP')*(spP*spX'), Npar ));
            D2L =  S1; % + S2;        

            Iinv = -(S1 - Sreg)^-1; 
            DI = Cxx*(repmat((1-2*P).*P.*(1-P),Npar,1).*X)';
            DL = DL + .5*DI'*Iinv(:) ;% + .5*DI2'*Iinv2(:); %Modified score function by Jeffreys prior (see Firth 1993)
            LL1 = sum(log(P*Yblocksum)) + .5*log(det(-D2L));                     

%             I1 = -D2L; I2 = -D2L; 
%             I1(:) = I1(:) + DI1*del1;% +  DI2*del1; %First order approximation to new I's
%             I2(:) = I2(:) + DI1*del2;% +  DI2*del2; %This might affect convergence.
            P =  exp( Theta2' * X )./  blocknorm(exp( Theta2' * X ) ,b); 
            dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
            spP = sparseblock(P,b);
            S1 = full( - X*dgP*X'  + spblocktrace( (spX*spP')*(spP*spX'), Npar ));
            D2L =  S1; % + S2;        

            Iinv = -(S1 - Sreg)^-1; 
            DI = Cxx*(repmat((1-2*P).*P.*(1-P),Npar,1).*X)';
            DL = DL + .5*DI'*Iinv(:) ;% + .5*DI2'*Iinv2(:); %Modified score function by Jeffreys prior (see Firth 1993)
            LL2 = sum(log(P*Yblocksum)) + .5*log(det(-D2L));                     

%             P =  exp( Theta2' * X )./  blocknorm(exp( Theta2' * X ) ,b); 
%             Iinv = -(S1 - Sreg)^-1; 
%             DI = Cxx*(repmat((1-2*P).*P.*(1-P),Npar,1).*X)';
%             DL = DL + .5*DI'*Iinv(:) ;% + .5*DI2'*Iinv2(:); %Modified score function by Jeffreys prior (see Firth 1993)
%             LL2  = sum(log(P*Yblocksum)) + .5*log(det(I2));
        end            


        if (LL > LL1 || isnan(LL1)) && ( LL > LL2 || isnan(LL2))
             damp = damp*lmnu;
 
        elseif LL1 > LL2 || isnan(LL2)
            
            Theta = Theta1;
            del = del1;
   
        else
            Theta = Theta2;
%             TH = TH2;
            del = del2;
            damp = damp./lmnu;
    
  
        end
        Slev = damp*eye(length(Theta));
%     elseif surrogate
%         
%         del1 = -D2L^-1  * (DL - Sreg*Theta);
%         Theta1 = Theta + del1;
%         P =  exp( Theta1' * X )./  blocknorm(exp( Theta1' * X ) ,b); 
%         LL1 = sum(log(P*Yblocksum));    
%         
%         if LL1 < LL
%              damp = damp*lmnu;
%              Slev = damp*eye(length(Theta));
%              D2Linv = (Hsurr - Slev)^-1;
%         else            
%             Theta = Theta1;
%             LL = LL1;
%             del = del1;
%         end
%         
    end

        
    
    dstep = sqrt( sum(del.^2)./sum(Theta.^2));
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    if ~isinf(dstep)
        fprintf('%5i: damping: %1.2e,  dstep = %1.2e',count,damp,dstep);
    else
        fprintf('%5i: damping: %1.2e,  dstep = %1.2e     ',count,damp,dstep);
    end
        
        count = count+1;

         dsteps(count) = dstep;
         
         llhist(count) = LL;

         if length(llhist)> 1& LL-llhist(end-1) <0, keyboard; end
         
             
         if count > maxcount % && std(llhist(end-50:end))./diff(llhist([1,end])) < 1e6 
            
             warning('Failed to converge after %i ieration. Estimate may be unbounded.',maxcount)
             pause(1)
             badcond = 1;
             break
         end
end
fprintf('\n');

%One more pass to Compute final value of everyting and observed Fisher information 
P =  exp( Theta' * X )./  blocknorm(exp( Theta' * X ) ,b); %Proabability of each outcome
LL = sum(log(P*Yblocksum));
Pnorm = blocknorm(P.*Y',b);
P2 = P.*Y'./Pnorm;
X2 = X(:,b2bl.*Y' > 1);
P2 = P2(:,b2bl.*Y' > 1);
dgP2 = sparseblockmex(P2,0:length(P2), 0:length(P2),length(P2)); %Equivalent to but faster than diag(P) or diag(sparse(P))
spP2 = sparseblock(P2,b2);
S2 = full( X2*dgP2*X2'  - spblocktrace( (spX2*spP2')*(spP2*spX2'), Npar ));
dgP = sparseblockmex(P,0:length(P), 0:length(P),length(P)); %Equivalent to but faster than diag(P) or diag(sparse(P))      
spP = sparseblock(P,b);
D2L = full( S2 - X*dgP*X'  + spblocktrace( (spX*spP')*(spP*spX'), Npar ));
if firth
    LL = LL + .5*log(det(-D2L));         
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

