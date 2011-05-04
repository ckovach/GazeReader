function varargout = setBinData(h,bindata)


%  setBinData(h,bindata)
%
% Associates the bindata with current data set
% and brings up the bin data window
%
% bindata is a structure returned by MAKEBINDATA
%
% See also makeBinData

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end



if nargin > 1
    setappdata(h,'binData',bindata);
end
    
if nargout > 0
    varargout{1} = getappdata(h,'binData');
end



handles = guidata(h);

GazeReader('binManagerMenu_Callback',h,[],handles);



emfun = getappdata(h,'EventManagerFunctions');
tmfun = getappdata(h,'trialManagerFunctions');

if ~isempty(emfun)
    emfun.updateAllTrials();
end

if ~isempty(tmfun)
    tmfun.updateAllDataSets();
end

if nargout > 1
    varargout{2} = getappdata(h,'binManagerFunctions');
end

