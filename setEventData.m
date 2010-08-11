function setEventData(h,eventdata)


%  setEventData(h,evtdata)
%
% Associates the evtdata with current data set
% and/or brings up the event data window
%
% evtdata is a structure returned by MAKEEVENTDATA
%
% See also MAKEEVENTDATA



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end



if nargin > 1
    setappdata(h,'expEventData',eventdata);
end

handles = guidata(h);

GazeReader('EventManager_Callback',h,[],handles); %Initialize events


