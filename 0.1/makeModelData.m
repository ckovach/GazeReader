
function modelDat = makeModelData(varargin)


label = '';
codeincr = 0;

if nargin > 0 && isstruct(varargin{1})
    modelDat = varargin{1};
    i = 2;
elseif nargin > 0 && isnumeric(varargin{1})
     i = 2;
else
    i = 1;
end

regressors = [];
Hreg = 0;
Lreg = 0;
while i <= length(varargin)
    
    switch lower(varargin{i})
        
        case 'label' %label
            label  = varargin{i+1};
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

            i = i+1; 
        case 'codeincr' %codeincrement
           codeincr = varargin{i+1};
            i = i+1; 
        case 'regressors' %codeincrement
           regressors = varargin{i+1};
            i = i+1; 
        case {'regularization','gaussreg'} %codeincrement
           Hreg = varargin{i+1};
            i = i+1; 
        case {'laplreg'}
           Hreg = varargin{i+1};
            i = i+1; 
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;    
end


% mdstruct = struct('label',label,'submodels',[],'parameterFit',[],'I',[],'badcond',[],...
%                                 'LL',[],'AIC',[],'BIC',[],'llrpval',[],'regressors',[],'contrasts',[],'code',codeincr+1);
mdstruct = struct('label',label,'fit',[],'regressors',regressors,'contrasts',[],'code',codeincr+1,'regLabels',[],'Hreg',Hreg,'Lreg',Lreg);

if nargin == 0 || ~isstruct(varargin{1})   
    modelDat = struct('models',[],'Y',[],'codeincr',0);
    modelDat.models = mdstruct;
else
    if isnumeric(modelDat.models)
        modelDat.models = mdstruct;
    else
        modelDat.models(end+1) = mdstruct;
    end
end

modelDat.codeincr = max([ modelDat.models.code]);
 

