function setTrialData(h,trialdata)


%  setTrialData(h,trialdata)
%
% Associates the trialdata with current data set
% and brings up the trial data window
%
% trialdata is a structure returned by MAKETRIALDATA
%
% See also makeTrialData



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end



if nargin > 1
    setappdata(h,'trialData',trialdata);
end

handles = guidata(h);

GazeReader('trialManagerMenu_Callback',h,[],handles);

