
function fit = modelFit(trialData ,R,varargin)

% function fit = modelFit(W,R,llrtests)
%Fits model. W is vector of outcomes. R is an array of regressor structures.
%Submodels are fit and LLR test carried out on all regressor groups in R if
%llrtests is true
%
% SEE MNLFIT

%C Kovach 2008

i=1;
Hreg = 0; %Default regularization
Firth = false;
fullOnly = false;
multassign =false;
LC = [];
include_null = 0;
binvolume = [];
llrtests = [];
inittheta = 0; %starting point for numerical maximization
% parallel = false; %use of parallel toolbox to be implemented
% maxcpu = 8; %Maximum workers used when using parallel computing
while i <= length(varargin)
    
    switch lower(varargin{i})
        case 'regularization'   %Specify gaussian regularization (ridge)
            Hreg = varargin{i+1};
            i = i+1;            
        case 'firth'    %Use Jeffrey's prior as described by Firth
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
        case 'multassign'     %Allow multiple assignment to bins
            multassign = varargin{i+1};
            i = i+1;
        case 'linearconstraint'     %A matrix of linear constraints on maximization
            LC = varargin{i+1};
            i = i+1;
        case 'inittheta'     
            inittheta = varargin{i+1};
            i = i+1;
        case 'binvolume'     %LOG bin volume
            binvolume = varargin{i+1};
            i = i+1;
        case 'parallel'     %use parallel computing toolbox if available
            parallel = varargin{i+1};
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
                                'LL',[],'LLR',[],'AIC',[],'BIC',[],'llrpval',[],'R',[],'regressors',[],'contrast',[],'W',[],'N',[],'firth',[],'Hreg',[],'multassign',[],'blockC',[],'discarded',[],'binvolume',[]);


                            
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
    W = sparse(zeros(sum(noptions),1));
    W(cumsum([0;noptions(1:end-1)]) + Wnum) = 1;
%     FB = Wnum;
    FB  = Wnum;
    FBr = sparseblock(ones(size(W)),noptions','transpose')*Wnum;
else
    W = double(trialData);
    if length(noptions) == 1;
            noptions = noptions*ones(1,length(W)./noptions);
    end
    FB = ones(size(noptions));
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


%Where W = 0, data are discarded
Rpooled.value = Rpooled.value(FBr~=0, any( Rpooled.value(FBr~=0,:)~=0,1) );
Worig = W;
W = W(FBr ~= 0);
Rpooled.noptions = Rpooled.noptions(FB ~= 0);

if length(binvolume) > 1
    binvolume = binvolume(FBr~=0);
end

nobs = sum(FB~=0);
% [parest,I,LL,badcond] = logitmodelsp(Rpooled,W, [] , Hreg);
% [parest,I,LL,badcond] = mnlfit(Rpooled,W,0, Hreg, true, [],Firth);
[parest,I,LL,badcond,lgm,max_iter] = mnlfit(Rpooled,W,'inittheta',inittheta,'regularization',Hreg, 'runiter',true, 'fix',[],'firth',Firth,...
                                'linearconstraint',LC,'include_null',include_null,'binvolume',binvolume,'checkdesign',true);

%The point estimate corresponding to fixed values and 0 otherwise
% [parest0,I0,LL0] = mnlfit(Rpooled,W,0, Hreg, false, [],Firth);
[parest0,I0,LL0] = mnlfit(Rpooled,W,'inittheta',inittheta, 'regularization',Hreg, 'runiter',false, 'fix',[],'firth',Firth,...
                                'linearconstraint',LC,'include_null',include_null,'binvolume',binvolume);

npar = length(parest) - rank(LC);
fit(1) = fitstruct;
fit(1).label = 'full';
fit(1).parest = parest;
fit(1).I = I;
fit(1).badcond = badcond;
fit(1).max_iterations_reached = max_iter;
fit(1).LL = LL;
fit(1).LLR = LL0-LL;
fit(1).llrpval = 1-chi2cdf(2*(LL - LL0),npar); %log likelihood ratio test against null hypothesis of no effect
fit(1).R = 1-LL./LL0; %Pseudo R-squared based on the ratio of the deviance of the model and the model with maximum entropy.
fit(1).npar = npar;
% fit(1).AIC = -2*LL + 2*npar; %Akaike Information Criterion
fit(1).AIC = -2*LL + 2*npar + 2*npar*(npar+1)./(length(Rpooled.noptions) - npar-1); %Akaike Information Criterion
                                                                                    % With second order correction
fit(1).BIC = -2*LL + npar*log(nobs); %Bayes-Schwartz information criterion
% fit(1).lBayes = npar./2*log(2*pi) - .5*log(det(I)) + LL; %Laplace's approximation to the log marginal likelihood. This can be used to compute approximate bayes factors

fit(1).contrast = [];%eye(sum([Rpooled.Npar]));
fit(1).regressors = [];
fit(1).W = Worig;
fit(1).blockC = blockC;
fit(1).firth = Firth;
fit(1).multassign= multassign;
fit(1).Hreg = Hreg;
fit(1).N = length(Rpooled.noptions);
fit(1).discarded = FBr == 0;
fit(1).binvolume = binvolume;

fullfit = fit;

regind = 1:length(R);
if length(R)>1 && ~fullOnly 
    
    fit((1:length(llrtests))+1) = fitstruct;
%     for i = 1:length(R)   

%     parfor (i = 1:length(R),maxcpu*parallel)  %Using parallel toobox
    regcodes = [R.code];
    for i = 1:length(llrtests)
        
        getreg = ~ismember(regcodes,llrtests{i});
        Rpooled = pool(R(getreg));
        
        Rpooled.value(isnan( Rpooled.value)) = 0;
        Rpooled.value = Rpooled.value(FBr~=0, any( Rpooled.value(FBr~=0,:)~=0,1) );
        Rpooled.noptions = Rpooled.noptions(FB ~= 0);

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
            LCsub = LCcontr'*LC;
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
        [parest,I,LL,badcond,lgm,max_iter] = mnlfit(Rpooled,W,'inittheta',stth, 'regularization',Hreg, 'runiter',true,...
                                                     'firth',Firth,'linearconstraint',LCsub,'binvolume',binvolume);

        npar = length(parest);
        fit(i+1).parest = parest;
        fit(i+1).I = I;
        fit(i+1).badcond = badcond;
        fit(i+1).max_iterations_reached = max_iter;
        fit(i+1).LL = LL;
        fit(i+1).LLR = LL-fullfit.LL;
        fit(i+1).R = fit(i+1).LLR./LL0; %Pseudo R-squared based on the ratio of the deviance of the model and the model 
                                       %with maximum entropy minus the same for the reduced model.

        fit(i+1).npar= npar;
%         fit(i+1).AIC = -2*LL + 2*npar;
        fit(i+1).AIC = -2*LL + 2*npar + 2*npar*(npar+1)./(length(Rpooled.noptions) - npar-1); %Akaike Information Criterion
                                                                                    % With second order correction
        fit(i+1).BIC = -2*LL + npar*log(nobs);
%         fit(i+1).lBayes = npar./2*log(2*pi) - .5*log(det(I)) + LL;
        fit(i+1).llrpval = 1-chi2cdf(2*(fullfit.LL - fit(i+1).LL),fullfit.npar-fit(i+1).npar); %log likelihood ratio test
        
        fit(i+1).regressors = [R(~getreg).code]; %Codes for the tested regressors
        
    end
end
