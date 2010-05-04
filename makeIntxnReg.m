
function regData = makeIntxnReg(varargin)

if isstruct(varargin{1}) && nargin > 1 && isstruct(varargin{2})
    regData = varargin{1};
    regData.regressors(end+1) = interaction(varargin{2:end},'codeincr',regData.codeincr);
    regData.codeincr =max([regData.regressors.code]);
else
    regData.regressors = interaction(varargin{1:end});
    regData.codeincr =max([regData.regressors.code]);
end
