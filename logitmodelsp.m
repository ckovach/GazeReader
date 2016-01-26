

function [Theta, I, LL, badcond] = logitmodelsp(R,Wsp,InitTheta, Sreg, runiter, fix, firth, surrogate )

%function [Theta, I] = logitmodelsp(SpX,Wsp,InitTheta, Sreg )
%
%     SpX is a block diagonal matrix. Each block, Bi is an b(i) x M matrix. M
%     is the number of regressors, which is fixed across trials. b(i) is the
%     number of possible outcomes on trial i, which varies from trial to trial.
%     Wsp is a sparse vector conatining 1 for each row in SpX which corresponds 
%     to an observed outcome.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

SpX = R.value;
b = R.noptions;
Nvar = R.Npar;

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
    InitTheta = zeros(Nvar,1) + eps;
elseif length(InitTheta ) < Nvar
    InitTheta(end+1:Nvar) = 0;
end

if nargin < 4 || isempty(Sreg)
    Sreg = 0; %Regularization term
elseif length(Sreg) == 1
    Sreg = eye(Nvar).*Sreg;
end

if isempty(SpX)
    runiter = false;
    InitTheta = [];
end

if LevMar
    damp = .1;
    Slev = damp*eye(length(InitTheta));

end

if transposeSpx 
    SpX = SpX';
end    

tol = 1e-6;


%Xs associated with the outcome
% XX = zeros(size(X,1),size(X,1),size(X,3),size(X,2));

nblocks = length(b);

SE = repmat(speye(Nvar), 1, nblocks );  % Used to compute the summation.

WX = SpX * Wsp;

% ncols =size(SpX,2);
% 
% [q,Ir,Jc] = unsparsifymex(SpX,'transpose');
% [Ir,Jc] = sparseij(SpX');
% clear q;

                                
% SEblk = repmat(speye(Nvar),nblocks,1);


Theta = sparse( InitTheta(:) );

dstep = Inf;
count = 0;
fprintf('%5i: damping: %1.2e,  dstep = %1.2e     ',0,lmnu,Inf);

% useloop = false;
% if Nvar^2*nblocks > maxMatSize
%     
%     useloop = true;
%     getblock(cumsum([1,b])) = 1;
%     getblock = cumsum(getblock(1:end-1));
%     getpar = kron(1:nblocks,ones(1,Nvar));
% end

del = Inf;


TH = repmat(diag(Theta),1,nblocks);
% LL = (  sum(TH,1) * SpX  - log(( exp( sum(TH,1) * SpX )*SEtrial)) ) * Wsp;
LL = (  sum(TH,1) * SpX  - log(blocknorm( exp( sum(TH,1) * SpX ),b) ) ) * Wsp;

while dstep > tol && runiter

%     THsp = sparseblockmex(repmat(theta,ncols,1),Ir,Jc); %Performance can be improved with a dedicated mex function here
    
    
    %probability  for each outcome
%     P =  exp( sum(TH,1) * SpX )./( exp( sum(TH,1) * SpX )*SEtrial); 
    P =  exp( sum(TH,1) * SpX )./blocknorm( exp( sum(TH,1) * SpX ), b); 
    
    % Derivative of the log likelihood
    
    wgtsum = SpX*P'; %Weighted sum 
    DL = SE*(WX - wgtsum) ;
    
    
    % Computing the second derivative
    

    wgtsumblk = sparseblock(full(wgtsum'),ones(1,nblocks)*Nvar);
    
    % terms for sum(  sum(  pi*Xi     ) * sum(  pi*Xi  )' ) and sum(  sum(  pi*Xi*Xi' )   )

    dgP = sparseblockmex(sqrt(full(P)),0:length(P), 0:length(P),length(P)); %faster than diag(P) or diag(sparse(P))

   

    D2L = (SE* wgtsumblk')*(wgtsumblk*SE')- (SE*SpX*dgP)*(SE*SpX*dgP)'; %This is the slow step
        
%     if firth
%         
%         H = dgP.*sqrt(1-dgP.^2)*spX*(SpX*dgP.^2.*(1-dgP.^2)*spX')^-1;
        
        
    % Newton's method
    
    if ~LevMar  
        del = -(D2L - Sreg)^-1 * (DL - Sreg*Theta);
        del(fix) = 0;

        Theta = Theta + del;
            % Computing the log likelihood
        TH = repmat(diag(Theta),1,nblocks);

    
    else
        
%         del1 = -(D2L - Sreg*lmnu)^-1 * DL;
        del1 = -(D2L - Slev - Sreg)^-1  * (DL - Sreg*Theta);
        del1(fix) = 0;
        del2 = -(D2L - Slev./lmnu - Sreg)^-1 * (DL - Sreg*Theta);
        del2(fix) = 0;
        Theta1 = Theta + del1;
        Theta2 = Theta + del2;
        

        TH1 = repmat(diag(Theta1),1,nblocks);
        TH2 = repmat(diag(Theta2),1,nblocks);

        LL1 = (  sum(TH1,1) * SpX  - log( blocknorm( exp( sum(TH1,1) * SpX ), b ) ) ) * Wsp;
        LL2 = (  sum(TH2,1) * SpX  - log( blocknorm( exp( sum(TH2,1) * SpX ), b ) ) ) * Wsp;


        if LL > LL1 && LL > LL2
             damp = damp*lmnu;
            %Sreg = Sreg*lmnu;
%             Theta = Theta;
        elseif LL1 > LL2
            
            Theta = Theta1;
            TH = TH1;
            del = del1;
            %Sreg = Sreg*lmnu;
            LL = LL1;
%             ext = ext1;
        else
            Theta = Theta2;
            TH = TH2;
            del = del2;
            damp = damp./lmnu;
            %Sreg = Sreg./lmnu;
            LL = LL2;
           
%             ext = ext2;
        end
%         LLs(end+1) = LL
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
%Log likelihood

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


