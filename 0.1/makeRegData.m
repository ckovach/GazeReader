
function regData = makeRegData(varargin)

if isempty(varargin(1))
    
    regData.regressors = makeregressor([],varargin{2:end});
    regData.regressors(1) = [];
    regData.condeincr = 0;
elseif isstruct(varargin{1})
    regData = varargin{1};
    regData.regressors(end+1) = makeregressor(varargin{2:end},'codeincr',regData.codeincr);
    regData.codeincr =max([regData.regressors.code]);
else
    regData.regressors = makeregressor(varargin{1:end});
    regData.codeincr =max([regData.regressors.code]);
end


