
function regData = makePolyReg(varargin)

if isstruct(varargin{1})
    regData = varargin{1};
    regData.regressors(end+1) = buildpolyreg(varargin{2:end},'codeincr',regData.codeincr);
    regData.codeincr =max([regData.regressors.code]);
else
    regData.regressors = buildpolyreg(varargin{1:end});
    regData.codeincr =max([regData.regressors.code]);
end
