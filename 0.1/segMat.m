
function varargout = segMat(ts,rg,fs)

%Create a matrix to segmant a signal sampled at fs round times ts over
%range rg.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 3
    fs = 1;
end

range = (rg(1):1/fs:rg(2))*fs;

T = round(repmat(ts(:)'*fs,length(range),1) + repmat(range(:),1,length(ts)));

T(:,any(T<1)) = [];

varargout{1} = T;

if nargout > 1
    varargout{2} = range/fs;
end

if nargout > 2
    varargout{3} = max(T(:));
end
