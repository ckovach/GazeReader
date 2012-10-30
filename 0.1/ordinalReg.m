
function varargout = ordinalReg(X,noptions,varargin)


% Create design matrix and response vector for an adjacent-categories
% ordinal logit (see Agresti 2010, chapter 3 & 4) for use with modelFit and
% mnlfit.
%
%  For adjacent-category logit, the probability for response i is
%   
%            Pi = Ai + beta*(i - bi)*X
%
%          where bi is some baseline category.
%
%  Rout = ordinalReg(Xin,noptions,'params',value);
%
%  Xin is an NxK matrix containing the indepenent variable as rows. 
%
%  Rout is a GazeReader regressor array with intercepts in the first
%  element, and the ordinal effect in the second. 
%
%  'noptions' is the total number of response categories which is either a
%  scalar if the number of options is the same for each row (as is
%  typically the case) or a column vector otherwise. 
%  
%  'noptions' can also be a  row vector. If so, then each column 
%  corresponds to a separate rating scale. Interaction terms between rating
%  scales are contained in the third element of Rout. In other words, each
%  outcome can be scored on multiple ordinal scales.
%  
%
%  The model is implemented as a single multinomial logit.
%  
%  By default the model assumes separate intercepts for each response
%  category. Alternatively the intercepts can be specified (as a function
%  of category number), by passing a function handle: 
%
%  Rout = ordinalReg(Xin,Yin,...,'index_fun', @(x) some_function(x),...);
%
% or intercepts can be disabled altogether
%
%          ,...,'intercepts', false,...);
%
% To specify any category as the baseline, use
%
%          ,...,'baseline', k,...); % Default k = 1;
%
%
%  [Rout, Rfull] = ordinalReg(...);
%
%  Returns the regressor for the full multinomial logit as a second
%  argument. For this the interaction with the input is modeled with a separate
%  set of parameters for each outcome category.
%
% Agresti A (2010) Analysis of ordinal categorical data. Hoboken, NJ: John Wiley & Sons Inc.
%
%  

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


baseline_index = 1; % This is the index of the options that serves as baseline.
cdi = 0;
label = 'covariates';
incl_intcpt = true;
indxfun = [];
poly_intcpt = 0;
i = 1;
%%%% Key words are described below %%%%
while i <= length(varargin)
    
    switch lower(varargin{i})
        case 'baseline'
            baseline_index = varargin{i+1};
            i = i+1;            
        case 'codeincr'
            cdi = varargin{i+1};
            i = i+1;            
        case 'label'
            label =  varargin{i+1};
            i = i+1;
        case 'intercepts'
            incl_intcpt =  varargin{i+1};
            i = i+1;
        case 'poly_intercept' %%% Use a polynomial of the given order
                              %%% for the intercepts rather than separate
                              %%% intercepts for each category
            poly_intcpt =  varargin{i+1};
            i = i+1;
            incl_intcpt = false;
        case 'index_fun'  %%% Alternative transformation on indices for the intercept term
            indxfun =  varargin{i+1};
            incl_intcpt = false;
            
            i = i+1;
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end

N = size(X,1);

if size(noptions,1) == 1
    noptions = repmat(noptions,N,1);
end


%%% Create a vector that expands each row of the design matrix to equal the
%%% number of options.

totnopt = prod(noptions,2);
endind = zeros(sum(totnopt),1);
trend =cumsum( totnopt(1:end));
endind(trend) = totnopt;
expand = cumsum([1;endind(1:end-1)~=0]);

%%% Create the GazeReader regressor for X
RX = makeregressor(X(expand,:),'noptions',totnopt,'label','X');


optn = [1:sum(totnopt)]' - cumsum([0;endind(1:end-1)]); %option number
if size(noptions,2) < 2
    indx = optn;
else
    
    %%% Note that noptions must be the same on every trial.
    %%% Variable noptions over trials is not implemented here.
    
    indx = repmat(crossn(noptions(1,:)),size(noptions,1),1);
    
end

if incl_intcpt
    %%% Create dummy variables for the options excluding the baseline
    RF = fact2reg(optn,'center',false,'ignore',optn == baseline_index,...
                  'noptions',totnopt,'label','intercepts','codeincr',cdi); cdi = cdi+1;
elseif poly_intcpt > 0 
        
    %%% Create polynomial basis functions over the indices of options
       RF = buildpolyreg(indx,poly_intcpt,'codeincr',cdi,'noptions',totnopt,'label','intercepts poly','codeincr',cdi); cdi = cdi+1;
       
elseif ~isempty(indxfun)
    for i = 1:size(indx,2)
       RF(i) = makeregressor(indxfun(indx(:,i)),'codeincr',cdi,'noptions',totnopt,'label',sprintf('intercepts %i',i),'codeincr',cdi); cdi = cdi+1;
    end
else
    RF = [];
end


%%% Compute the interaction with option number. This is equivalent to the 
%%% assumption of the adjacent categories logit:
%%%
%%%      log(p_i/p_{i-1}) = (a_i-a_{i-1}) + b*x
%%% 
%%%  It follows
%%%
%%%      log(p_i/p_0) = a_i + i*b*x
%%%
%%%  where p_0 is the probability for the baseline category.



%%% Regressor for the option number (used
Rnopt = makeregressor(indx-baseline_index,'noptions',totnopt,'label','order');

%%% This allows for multiple simultaneous ordinals
if size(indx,2) > 1
    Rnopt = split(Rnopt);
end
%%% INteraction term
Rintxn = interaction(Rnopt,RX,'codeincr',cdi);cdi = cdi+1;

if length(Rnopt) > 1
    Rintint = interaction(Rnopt,'codeincr',cdi);cdi = cdi+1;
    Rintintx = interaction(Rintint,RX,'codeincr',cdi);cdi = cdi+1;
else
    Rintintx = [];
    Rintxn = [];
end
%%% Array of output regressors
varargout{1} = [RF,Rintxn, Rintintx]; 


%%% If more than one output, then make the second set of regressors for the
%%% non-ordinal full multinomial logit.

if nargout > 1    
       varargout{2} = [RF, interaction(RF,RX,'label','full multinomial','codeincr',cdi)]; cdi = cdi+1;
             
end
           



function tab = crossn(N,full)

% tab = crossn(N)
%
% Creates a matrix whose rows contain all combinations of integers
% 1:N(1),1:N(2),....
%
% For example, if N is a 2 X 1 vector of the lengths of 2 vectors, x1 and x2,
% then crossn creates a  N(1)*N(2) X 2 matrix which contains every pair of
% indices into x1 and x2.
%
%
% tab = crossn(N,0)
%
% Returns only combinations that are unique under permutation.


% C. Kovach 2010


if nargin < 2 || isempty(full)
    full = true;
end


if length(N) == 1    
    tab = (1:N)';
    return
end
    
cn = crossn(N(2:end),full);
    
tab = cat(2,kron((1:N(1))',ones(size(cn,1),1)), repmat(cn,N(1),1));    
    
if ~full %Discard permutations of the same vector
    
%    [sqi,sqj]  = meshgrid(1:N(1),1:size(cn,1)); 
%    
%    keep = sqi <= sqj | sqi > size(cn,1);
   
   cols = 1:size(tab,2);
   
   [mn,mni] = min(max(tab));        

   keep = all(repmat(tab(:,mni),1,length(cols)-1) - tab(:,setdiff(cols,mni)) <=0,2);  
   
   tab = tab(keep,:);

end

