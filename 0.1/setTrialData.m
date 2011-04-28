function trialdata = setTrialData(h,trialdata)


%  setTrialData(h,trialdata)
%
% Associates the trialdata with current data set
% and brings up the trial data window
%
% trialdata is a structure returned by MAKETRIALDATA
%
% See also makeTrialData

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
    setappdata(h,'trialData',trialdata);
elseif nargout > 0
    trialdata = getappdata(h,'trialData');
end
handles = guidata(h);

GazeReader('Trial_Manager_Callback',h,[],handles);

tmfun = getappdata(h,'trialManagerFunctions');

if ~isempty(tmfun)
      tmfun.updateAllDataSets();
 %   tmfun.update();
end


