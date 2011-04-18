
function fit = modelFit(trialData ,R,varargin)

% function fit = modelFit(Y,R,option,'value')
%
% Function to fit multinomial generalized linear models and compute
% likelihood ratio tests.
%
% R is an array of regressor group structures as returned by makeregressor. 
%
% Y is the vector of observations. 
%
% Submodels are fit and likelihood ratio tests are carried out on all 
% regressor groups in R by default.
%
% fit is a 1 x K structure where fit(1) contains data on the model fit
% for the full model, and fit(i) contains data for the (i-1)th likelihood
% ratio test.
%  
%    .label: label for the test
%
%    .parest: model parameter estimates for the full model in fit(1), and
%            for the model excluding the (i-1)th group of regressors otherwise. 
%
%    .npar:  number of parameters in the fitted model.
%
%    .I   :   Observed Fisher information matrix.
%
%    .badcond: 1 if the Fisher information matrix is singular, 0 otherwise.
%
%    .max_iterations_reached: 1 if fitting terminated after the maximum
%                number of iterations was reached.
%     
%    .LL:     log likelihood for the parameter estimates.
% 
%    .LLR:    log likelihood ratio with respect to the full model. For
%           fit(1).LLR, LLR is the likelihood ratio with respect to the maximum
%           entropy model (having uniform probability). Note that deviance
%           is given by -2*fit(i).LLR.
%           
%    .AIC:  Bias corrected Akaike's information criterion.
% 
%    .BIC:  Bayes information criterion.
% 
%    .llrpval: p-value from the likelihood ratio test, which refers
%            deviance (-2*LLR) to a chi-square distribution with degrees of
%            freedom equal to the difference in the number of parameters
%            between the full and reduced model.
% 
%    .R:    pseudo R-squared statistic.
%
%    .regressors: codes for the regressors included in the LLR test.
%
%    .contrast: Contrast matrix to extract the regressors included in the LLR test. 
%
%    .W: Vector of observation
%
%    .N: number of data points.
%
%    .obsfreq: frequency of each observation. Model fitting is more
%             efficient when observations that are repeated with respect to design 
%             matrix and outcome are represented by a single row in the
%             design matrix with the number of occurrences in obsfreq.
%   
%    .firth: 1 if regularization with Jeffrey's prior is used. This currently works only with
%           dichotomous outcomes.: 
% 
%    .Hreg:  inverse covariance matrix of a Gaussian prior distribution on
%       the parameters.
%
%    .Multassign: 1 if multiple assignment to bins is treated as a mixture
%            model.
%    
%    .discarded:  Indicates discarded data points
%    
%    .binvolume:  the volume (area) of each bin.
%  
%    .singular_design_matrix: True if design matrix is singular.
%
%  Options: 
%          'regularization'    Specify gaussian prior covariance matrix
%          'firth'             Use Jeffrey's prior as described by Firth -- this currently works only for dichotomous logistic regression with one 
%          'include_null'      Include a constant term for a "null" option not explicitly represented 
%                                   in the outcome vector and regressor matrix.
%          'fullonly'          Only fits the full model (no LLR statistic on submodels
%          'llrtests'          a cell array of groups of regressors on which to conduct likelihood ratio tests.
%          'multassign'        Allow multiple assignment to bins ( mixture models )
%          'linearconstraint'  A matrix of linear constraints on maximization
%          'inittheta'         Starting value of the parameter estimate.
%          'H0theta'           value of the parameters for the null hypothesis in the global test (defaults to 0, implying maximum entropy)   
%          'binvolume'         LOG bin volume
%          'obsfreq'           Vector of observation frequncies for each
%                              observation in the input.
%           'discard'     	   True for every data point which should eb
%           'fix'              1 x npar vector which is non-zer0 for every regressor whos parameter estimate 
%                              should remain fixed at the specified value.
%          {'diagsonly','show_progress','maxiter'}     %Other options passed to mnlfit
%
%
% SEE MAKEREGRESSOR, MNLFIT

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C Kovach 2008

i=1;
Hreg = 0; %Default regularization
Firth = false;
fullOnly = false;
multassign =true;
LC = [];
include_null = 0;
binvolume = [];
llrtests = [];
obsfreq = 1;
discard = 0;
inittheta = 0; %starting point for numerical maximization
collapserows = false;  %Improve efficiency by calling collapseX to collapse over like rows
% parallel = false; %use of parallel toolbox to be implemented
% maxcpu = 8; %Maximum workers used when using parallel computing


H0theta = 0;

fix = [R.fixed]; %specified parameters remain fixed
mnlfitopts = {};

while i <= length(varargin)
    
    switch lower(varargin{i})
        case 'regularization'   %Specify gaussian regularization (ridge)
            Hreg = varargin{i+1};
            i = i+1;            
        case 'firth'    %Use Jeffrey's prior as described by Firth -- this currently works only for dichotomous logistic regression with one 
            Firth = varargin{i+1};
            i = i+1;
       case 'include_null'    %Include a constant term for a "null" option not explicitly represented 
                              % in the outcome vector and regressor matrix.
            include_null  = varargin{i+1} ;
            i = i+1;
        case 'fullonly'     %Only fits the full model (no LLR statistic on submodels
            fullOnly = varargin{i+1};
            i = i+1;
        case 'llrtests'     %groups of regressors on which to conduct log-likelihood ratio tests.
            llrtests = varargin{i+1};
            i = i+1;
            fullOnly = false;
        case 'multassign'     %Allow multiple assignment to bins ( mixture models )
            multassign = varargin{i+1};
            i = i+1;
        case 'linearconstraint'     %A matrix of linear constraints on maximization
            LC = varargin{i+1};
            i = i+1;
        case 'inittheta'     
            inittheta = varargin{i+1};
            i = i+1;
        case 'h0theta'  %value of the parameters for the null hypothesis in the global test (defaults to 0)   
            H0theta = varargin{i+1};
            i = i+1;
        case 'binvolume'     %LOG bin volume
            binvolume = varargin{i+1};
            i = i+1;
       case 'obsfreq'     %LOG bin volume
            obsfreq = varargin{i+1};
            collapserows = false; %Assume already collapsed.
            i = i+1;
       case 'discard'     %LOG bin volume
            discard= varargin{i+1};
            i = i+1;
        case 'parallel'     %use parallel computing toolbox if available - NOT CURRENTLY FUNCTIONAL
%             parallel = varargin{i+1};
            i = i+1;
        case 'fix'    
            fix = varargin{i+1};
            i = i+1;
        case {'diagsonly','show_progress','maxiter'}     %Other options passed to mnlfit
            mnlfitopts(end+(1:2)) = varargin(i+(0:1)); 
            i = i+1;
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end



noptions = R(1).noptions;

if any(diff(cat(2,R.noptions),[],2) ~= 0)
    error('Structure of trials is not the same for all regressors.')
end

% fitstruct = struct('label','','parest',[],'npar',[],'I',[],'badcond',[],...
%                                 'LL',[],'AIC',[],'BIC',[],'lBayes',[],'llrpval',[],'regressors',[],'contrast',[],'W',[],'firth',[],'Hreg',[],'multassign',[],'blockC',[]);
fitstruct = struct('label','','parest',[],'npar',[],'I',[],'badcond',[],'max_iterations_reached',[],...
                                'LL',[],'LLR',[],'Dev',[],'deltaDev',[],'df',[],'llrpval',[],'AIC',[],'BIC',[],'R',[],'regressors',[],...
                                'contrast',[],'W',[],'N',[],'obsfreq',[],'firth',[],'Hreg',[],'multassign',[],'blockC',[],...
                                'discarded',[],'binvolume',[],'singular_design_matrix',[]);

%%% Create the response vector
                            
if isstruct(trialData)
    
    % get binvolume data
    if isempty(binvolume) && isfield(trialData.trials,'binareas')
         bv = cellfun(@(X,a) repmat(X,a,1)',{trialData.trials.binareas},{trialData.trials.nfix},'uniformoutput',0);
        bv = cellfun(@(X) X(:),bv,'uniformoutput',0);
        binvolume = log(cat(1,bv{:}));
    end
    if isempty(binvolume), binvolume = 0; end
    

    %     FB = []; 
    FBr= []; % FxNum=[];
    nfixs = [trialData.trials.nfix];
    nbins = [trialData.trials.nbin];
    W = spalloc(sum(nbins.*nfixs),1,sum(nfixs));
    currind = 0;
    for k = 1:length(trialData.trials);
        nfix = trialData.trials(k).nfix;
         nbin = trialData.trials(k).nbin;    
        fb = trialData.trials(k).fixbin;
        
        if nfix == 0
            continue
        end
        if ~multassign 
            w = spalloc(nfix*nbin,1,nfix);
            c = (0:length(fb))';
            w( fb(fb ~=0) +c(fb~=0)*nbin) = 1;
%             W = [W;w];
            W(currind+(1:length(w))) = w;
            currind = currind + length(w);
        else
            fm = trialData.trials(k).fixmat';
            W(currind+(1:nfix*nbin)) = fm(:);
            currind = currind + nfix*nbin;
%             W = [W;fm(:)];
        end
%          FB = cat(1,FB, fb*ones(nfix*nbin,1) );
        FBr = cat(1,FBr, kron(fb,ones(nbin,1) ));
    %    FxNum = cat(1,FxNum,kron((1:nfix)',ones(nbin,1)) + sum([trialData.trials(1:k-1).nfix]));
    end
    FB = cat(1,trialData.trials.fixbin);

elseif isnumeric(trialData) && ~issparse(trialData)
    
    Wnum = trialData;
    if length(noptions) == 1;
            noptions = noptions*ones(1,length(Wnum));
    end
%     W = spalloc(sum(noptions),1,length(noptions));
    W(sum(noptions)) = 0;
    W(cumsum([0;noptions(1:end-1)]) + Wnum) = 1;
    W = sparse(W');
%     FB = Wnum;
    FB  = Wnum;
    FBr = sparseblock(ones(size(W)),noptions','transpose')*Wnum;
else
    W = double(trialData);
    if length(noptions) == 1;
            noptions = noptions*ones(1,length(W)./noptions);
    end
    
%     FB = ones(size(noptions));
    
    %Get the bin index for each event.
    numvec = zeros(size(W));
    numvec(cumsum(noptions(1:end))) = noptions;
    numvec = cumsum(numvec);
%    sbl = sparseblock((1:sum(noptions))-numvec',noptions); %modified   2/1/2010
    sbl = sparseblock((1:sum(noptions))-numvec'+1,noptions);
    FB = sbl*W;
%     FBr = ones(size(W));
    FBr = sparseblock(ones(size(W)),noptions','transpose')*FB;
end

if sum(noptions) ~= length(W)
    error('Sum of noptions does not match the size of the observation vector')
end


if isnumeric(R)
    R = makeregressor(R,'noptions',noptions);
end

if isempty(llrtests) && ~fullOnly
    llrtests = num2cell([R.code]);
end

Rpooled = pool(R);

Npars = [R.Npar];

blockC = full(sparseblock(ones(1,Rpooled.Npar),Npars)');
Rpooled.value(isnan(Rpooled.value)) = eps;


if any(fix)
    fixed = find(fix);
    LC = zeros(sum(Npars),1);
    for i = 1:length(fixed)
        LC(fixed(i),i) = 1;
    end
    if any(inittheta~=0)
        LC(end+1,:) = inittheta(fixed);
    end
end

%Where W = 0, data are discarded
% Rpooled.value = Rpooled.value(FBr~=0, any( Rpooled.value(FBr~=0,:)~=0,1) );
Rpooled.value = Rpooled.value(:, any( Rpooled.value(FBr~=0,:)~=0,1) );
% Rpooled.noptions = Rpooled.noptions(FB ~= 0);
Worig = W;
% W = W(FBr ~= 0);

discard = discard | FB==0;

if length(binvolume) > 1
    binvolume = binvolume(FBr~=0);
end

nobs = sum(FB~=0);

if collapserows
    %%% collapse over like rows if collapserows is true
    %%% This can make fitting more efficient by converting binary outcomes 
    %%% to counts
  
    
    [colindx, obsfreq, expindx, nopt] = collapseX(cat(2,W,Rpooled.value),Rpooled.noptions,'discardtrials',discard);
    
    Rpooled.value = Rpooled.value(colindx,:);
    Rpooled.noptions = nopt;
    
        
    Win = W(colindx);
    
else 
%     obsfreq = ones(size(Rpooled.noptions));
    Win = W;
end
        
% [parest,I,LL,badcond] = logitmodelsp(Rpooled,W, [] , Hreg);
% [parest,I,LL,badcond] = mnlfit(Rpooled,W,0, Hreg, true, [],Firth);
[parest,I,LL,badcond,lgm,max_iter,design_singular,Nobs] = mnlfit(Rpooled,Win,'inittheta',inittheta,'regularization',Hreg, 'runiter',true, 'fix',[],'firth',Firth,...
                                'linearconstraint',LC,'include_null',include_null,'binvolume',binvolume,'checkdesign',true,'discard',discard,'obsfreq',obsfreq,mnlfitopts{:});

%The point estimate corresponding to fixed values and 0 otherwise
% [parest0,I0,LL0] = mnlfit(Rpooled,W,0, Hreg, false, [],Firth);
[parest0,I0,LL0] = mnlfit(Rpooled,Win,'inittheta',H0theta, 'regularization',Hreg, 'runiter',false, 'fix',[],'firth',Firth,...
                                'linearconstraint',LC,'include_null',include_null,'binvolume',binvolume,'discard',discard,'obsfreq',obsfreq,mnlfitopts{:});

npar = length(parest) - rank(LC);
fit(1) = fitstruct;
fit(1).label = 'full';
fit(1).parest = parest;
fit(1).I = I;
fit(1).badcond = badcond;
fit(1).singular_design_matrix = design_singular;
fit(1).max_iterations_reached = max_iter;
fit(1).LL = LL;
fit(1).Dev = -2*LL;
fit(1).LLR = LL0-LL;
fit(1).df = npar;
fit(1).deltaDev = -2*(LL0-LL);
fit(1).llrpval = 1-chi2cdf(2*(LL - LL0),npar); % log likelihood ratio test against null hypothesis of no effect
fit(1).R = 1-LL./LL0; % Pseudo R-squared based on the ratio of the deviance of the model and the model with maximum entropy.
fit(1).npar = npar;
% fit(1).AIC = -2*LL + 2*npar; %Akaike Information Criterion
fit(1).AIC = -2*LL + 2*npar + 2*npar*(npar+1)./(length(Rpooled.noptions) - npar-1); % Akaike Information Criterion
                                                                                    % With second order correction
fit(1).BIC = -2*LL + npar*log(nobs); % Bayes-Schwartz information criterion
% fit(1).lBayes = npar./2*log(2*pi) - .5*log(det(I)) + LL; % sLaplace's approximation to the log marginal likelihood. This can be used to compute approximate bayes factors

fit(1).contrast = [];%eye(sum([Rpooled.Npar]));
fit(1).regressors = [];
fit(1).W = Worig;
fit(1).blockC = blockC;
fit(1).firth = Firth;
fit(1).multassign= multassign;
fit(1).Hreg = Hreg;
% fit(1).N = length(Rpooled.noptions);
fit(1).N = full(Nobs);
fit(1).obsfreq = obsfreq;

fit(1).discarded = FBr == 0;
fit(1).binvolume = binvolume;

fullfit = fit;

% regind = 1:length(R);

if length(R)>1 && ~fullOnly 
    
    fit((1:length(llrtests))+1) = fitstruct;
%     for i = 1:length(R)   

%     parfor (i = 1:length(R),maxcpu*parallel)  %Using parallel toobox
    regcodes = [R.code];
    for i = 1:length(llrtests)
        
        getreg = ~ismember(regcodes,llrtests{i});
        Rpooled = pool(R(getreg));
        

        if collapserows
            %%% collapse over like rows if collapserows is true
            %%% This can make fitting more efficient by converting binary outcomes 
            %%% to counts


            [colindx, obsfreq, expindx, nopt] = collapseX(cat(2,W,Rpooled.value,discard),Rpooled.noptions);

            Rpooled.value = Rpooled.value(colindx,:);
            Rpooled.noptions = nopt;


            Win = W(colindx);

        else 
            Win = W;
        end
        
%         Rpooled.value(isnan( Rpooled.value)) = 0;
%         Rpooled.value = Rpooled.value(FBr~=0, any( Rpooled.value(FBr~=0,:)~=0,1) );
%         Rpooled.noptions = Rpooled.noptions(FB ~= 0);

%         fit(i+1).contrast = zeros(sum([R.Npar]),R(i).Npar);
%         fit(i+1).contrast( sum([R(1:i-1).Npar]) +(1:R(i).Npar),: ) = eye(R(i).Npar);
       
        rstartind = [0,cumsum([R(1:end-1).Npar])]+1;
        rendind = [0,cumsum([R(1:end-1).Npar])]+[R(1:end).Npar];
        Q = zeros(1,sum([R.Npar])+1);
        Q(rstartind(~getreg)) = 1;
        Q(rendind(~getreg)+1) = Q(rendind(~getreg)+1) -1;
        C = diag(cumsum(Q(1:sum([R.Npar]))));
        fit(i+1).contrast = C(:,sum(C)>0);
        
        if ~isempty(LC)
            invC = eye(size(C))-C;
            LCcontr = invC(:,sum(invC)>0);
            LCcontr(sum(Npars)+1,:) = 0;
            LCsub = LCcontr'*LC;
            LCsub(end+1,:) = LC(end,:);
            LCsub(:,sum(LCsub(1:end-1,:))==0)= [];
        else
            LCsub = [];
        end
        
        if size(LC,1)>size(C,2)
            LCcontr(end+1,end+1) = 1;
        end
            stth = 0;
%         fit(i+1).label = sprintf('%s',R(i).label);
        if sum(~getreg)==1
          fit(i+1).label = sprintf('%s',R(~getreg).label);
        else
            fit(i+1).label = sprintf(' | %s ',R(~getreg).label);
        end
        fprintf('\nFitting submodel %i: %s\n',i,fit(i+1).label );
        [parest,I,LL,badcond,lgm,max_iter] = mnlfit(Rpooled,Win,'inittheta',stth, 'regularization',Hreg, 'runiter',true,...
                                                     'firth',Firth,'linearconstraint',LCsub,'binvolume',binvolume,'discard',discard,'obsfreq',obsfreq,mnlfitopts{:});

        npar = length(parest) - rank(LCsub);
        fit(i+1).parest = parest; %#ok<*AGROW>
        fit(i+1).I = I;
        fit(i+1).badcond = badcond;
        fit(i+1).max_iterations_reached = max_iter;
        fit(i+1).LL = LL;
        fit(i+1).Dev = -2*LL;
        fit(i+1).LLR = LL-fullfit.LL;
        fit(i+1).deltaDev = -2*(LL-fullfit.LL);
        fit(i+1).R = fit(i+1).LLR./LL0; %Pseudo R-squared based on the ratio of the deviance of the model and the model 
                                       %with maximum entropy minus the same for the reduced model.
                                       
        fit(i+1).df = fullfit.npar-npar;
        
        fit(i+1).npar= npar;
%         fit(i+1).AIC = -2*LL + 2*npar;
        fit(i+1).AIC = -2*LL + 2*npar + 2*npar*(npar+1)./(length(Rpooled.noptions) - npar-1); %Akaike Information Criterion
                                                                                              % With second order correction
        fit(i+1).BIC = -2*LL + npar*log(nobs);
%         fit(i+1).lBayes = npar./2*log(2*pi) - .5*log(det(I)) + LL;
        fit(i+1).llrpval = 1-chi2cdf(2*(fullfit.LL - fit(i+1).LL),fullfit.npar-fit(i+1).npar); %log likelihood ratio significance test
        
        fit(i+1).regressors = [R(~getreg).code]; %Codes for the tested regressors
        
    end
end
