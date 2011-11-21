
function Rout = ordinalReg(X,noptions,varargin)

% Create design matrix and response vector for an adjacent-categories
% ordinal logit (see Agresti 2010, chapter 3 & 4) for use with mnlfit.
%
%  For adjacent-category logit, the probability for response i is
%   
%            Pi = Ai + beta*(i - bi)*X
%
%          where bi is some baseline category.
%
%  [Xout,Yout] = ordinalReg(Xin,Yin,noptions);
%
%  Xin is an NxK matrix containing the indepenent variable as rows. 
%
%  'noptions' is the total number of response categories which is either a
%  scalar if the number of options is the same for each row (as is
%  typically the case) or a vector otherwise. 
%
%  The model is implemented as a single multinomial logit.
%  
%  By default the model assumes separate intercepts for each response
%  category. Alternatively the intercepts can be specified (as a function
%  of category number), by passing a function handle: 
%
%  [Xout,Yout] = ordinalReg(Xin,Yin,...,'index_fun', @(x) some_function(x),...);
%
% or intercepts can be disabled altogether
%
%  [Xout,Yout] = ordinalReg(Xin,Yin,noptions,...,'intercepts', false,...);
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
index_fun = [];

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

if length(noptions) == 1
    noptions = repmat(noptions,N,1);
end


%%% Create a vector that expands each row of the design matrix to equal the
%%% number of options.

endind = zeros(sum(noptions),1);
trend =cumsum( noptions(1:end));
endind(trend) = noptions;
expand = cumsum([1;endind(1:end-1)~=0]);

%%% Create the GazeReader regressor for X
RX = makeregressor(X(expand,:),'noptions',noptions,'label','X');


optn = [1:sum(noptions)]' - cumsum([0;endind(1:end-1)]); %option number

if incl_intcpt
    %%% Create dummy variables for the options excluding the baseline
    RF = fact2reg(optn,'center',false,'ignore',optn == baseline_index,...
                  'noptions',noptions,'label','intercepts','codeincr',cdi); cdi = cdi+1;
elseif ~isempty(indxfun)
    RF = makeregressor(indxfun(optn),'codeincr',cdi,'noptions',noptions,'label','intercepts','codeincr',cdi); cdi = cdi+1;
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
Rnopt = makeregressor(optn-1,'noptions',noptions,'label','order');

%%% INteraction term
Rintxn = interaction(Rnopt,RX,'codeincr',cdi);cdi = cdi+1;

%%% Array of output regressors
Rout = [RF,Rintxn]; 



