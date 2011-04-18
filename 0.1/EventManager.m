function varargout = EventManager(varargin)
% EVENTMANAGER M-file for EventManager.fig
%      EVENTMANAGER, by itself, creates a new EVENTMANAGER or raises the existing
%      singleton*.
%
%      H = EVENTMANAGER returns the handle to a new EVENTMANAGER or the handle to
%      the existing singleton*.
%
%      EVENTMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVENTMANAGER.M with the given input arguments.
%
%      EVENTMANAGER('Property','Value',...) creates a new EVENTMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EventManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EventManager_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% Edit the above text to modify the response to help EventManager

% Last Modified by GUIDE v2.5 02-Jan-2008 14:40:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EventManager_OpeningFcn, ...
                   'gui_OutputFcn',  @EventManager_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EventManager is made visible.
function EventManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EventManager (see VARARGIN)

% Choose default command line output for EventManager
handles.output = hObject;

parent = varargin{1};
setappdata(handles.figure1,'parent',parent);
ch = getappdata(parent,'children'); ch(end+1) = handles.figure1;
setappdata(parent,'children',ch);

setappdata(parent,'EventManager',handles.figure1);

expEventData = getappdata(parent,'expEventData');
% currentDataSet = getappdata(parent,'CurrentDataSet');

evtfuns.update = @() Update(handles.figure1,[],handles);
evtfuns.updateTrials = @() UpdateTrialData(handles.figure1,[],handles);
setappdata(parent,'EventManagerFunctions',evtfuns);

if isempty(expEventData)
   errordlg('There appear to be no events. Assign trial onset and end in EventManager first.');
end

% if ~isfield(expEventData(currentDataSet),'events') || isempty(expEventData(currentDataSet).events)
% 
%     xdatTs = [expEventData(currentDataSet).xdat.startT];
%     xdatcodenum = [expEventData(currentDataSet).xdat.id];
%     type = zeros(1,length(xdatTs));
%     xdatlabel = {expEventgetappdata(h,'EventManager')Data(currentDataSet).xdat.code};
%     
%     
%     newevents = makeEventData(expEventData(currentDataSet),'time',xdatTs,'type',type,'xdatcode',xdatcodenum,'label',xdatlabel);
%     
%     expEventData(currentDataSet).events = newevents.events;    
%     expEventData(currentDataSet).codeincr = length(xdatTs);
%     
% end    
    
setappdata(parent,'expEventData',expEventData)



    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EventManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);

Update(hObject, eventdata, handles)
% --- Outputs from this function are returned to the command line.
function varargout = EventManager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in eventList.
function eventList_Callback(hObject, eventdata, handles)
% hObject    handle to eventList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns eventList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eventList

parent = getappdata(handles.figure1,'parent');

trmfun = getappdata(parent,'trialManagerFunctions');

if ~strcmp(get(handles.figure1,'selectionType'),'open')
    Update(hObject, eventdata, handles);
    if ~isempty(trmfun)
        trmfun.update(parent);
    end
    return
end

expEventData = getappdata(parent,'expEventData');
currentDataSet = getappdata(parent,'CurrentDataSet');
newEvtType = get(handles.eventTypeMenu,'value')-1;

selected = get(handles.eventList,'value')-1;

evtcodes = [expEventData(currentDataSet).events.xdatcode];
% evttypes= [expEventData(currentDataSet).events.type];
evttimes= double([expEventData(currentDataSet).events.time]);
evtlabels = {expEventData(currentDataSet).events.label};


if get(handles.markAllCodesCheck,'value');    
    selected = find(evtcodes== evtcodes(selected));
end

offset = str2num(get(handles.offset,'string'));

% if offset == 0  
%     c = mat2cell(newEvtType*ones(size(selected)));    
%     [expEventData.events.type] = c{:};        
% else
% 

% typestr = get(handles.eventTypeMenu,'string');

newevttimes = evttimes(selected) + offset;
newevttypes = newEvtType*ones(size(selected));
newevtlabels = evtlabels(selected); 
% newevtlabels = typestr(newevttypes+1);
newevtcodes = evtcodes(selected);
if get(handles.duplicateCheckBox,'value')
    expEventData(currentDataSet) = makeEventData(expEventData(currentDataSet),'type',newevttypes,'time',...
                newevttimes,'label',newevtlabels,'xdatcode',newevtcodes);
else
    idcodes = {expEventData(currentDataSet).events(selected).code};
    newevents =  makeEventData('type',newevttypes,'time',...
                newevttimes,'label',newevtlabels,'xdatcode',newevtcodes,'code',idcodes);
    expEventData(currentDataSet).events(selected) = newevents.events;
end
setappdata(parent,'expEventData',expEventData)
% end
if ismember(newEvtType,[1 2])
    UpdateTrialData(hObject,eventdata,handles)
end
Update(hObject, eventdata, handles);

set(handles.offset,'string','0')
% --- Executes during object creation, after setting all properties.
function eventList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double


% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function Update(hObject,eventData,handles)


% UpdateTrialData(hObject,eventData,handles)

parent = getappdata(handles.figure1,'parent');

expEventData = getappdata(parent,'expEventData');


currentDataSet = getappdata(parent,'CurrentDataSet');

% trialData = getappdata(parent,'trialData');

% if isempty(currentDataSet) || currentDataSet==0
%     currentDataSet = 1;
%     setappdata(parent,'CurrentDataSet',currentDataSet )
% end
if isempty(currentDataSet) || currentDataSet==0
    fprintf('No data set selected...')
    return
end


% expEventData = expEventData(currentDataSet);

if isempty(expEventData) || isempty([expEventData.events])
    error('There appears to be no loaded data set or no events.')
end

%String for the event Data list box
% evttypenum = cellfun(@num2str,{expEventData.xdat.id},'uniformoutput',0);

evtxdat = cellfun(@num2str,{expEventData(currentDataSet).events.xdatcode},'uniformoutput',0);
cstrlen = cellfun(@length,evtxdat);
spaces1 = cellfun(@(n) repmat(' ',1,n),num2cell(10-2*cstrlen),'UniformOutput',false);

typestr = get(handles.eventTypeMenu,'string');

typestr{strcmp(typestr,'Nothing')} = '';

evttypelbl = typestr([expEventData(currentDataSet).events.type]+1)';
cstrlen = cellfun(@length,evttypelbl);
spaces2 = cellfun(@(n) repmat(' ',1,n),num2cell(30-2*cstrlen),'UniformOutput',false);

evtlabels = {expEventData(currentDataSet).events.label};
cstrlen = cellfun(@length,evtlabels);
spaces3 = cellfun(@(n) repmat(' ',1,n),num2cell(12-2*cstrlen),'UniformOutput',false);

evt = [expEventData(currentDataSet).events.time];
evttimes = cellfun(@num2str,num2cell(evt),'uniformoutput',0);
tstrlen = cellfun(@length,evttimes);
spaces4 = cellfun(@(n) repmat(' ',1,n),num2cell(15-2*tstrlen),'UniformOutput',false);

number = cellfun(@num2str,num2cell(1:length(evt)),'uniformoutput',0);
tstrlen = cellfun(@length,number);
spaces5 = cellfun(@(n) repmat(' ',1,n),num2cell(10-2*tstrlen),'UniformOutput',false);

liststr = strcat(number,spaces5,evttimes,spaces4,evtxdat,spaces1,evttypelbl,spaces2,evtlabels,spaces3);

liststr = [{'No.  Time    code#     event type          code  '},liststr];

set(handles.eventList,'string',liststr)





%-------------------------------------------

function UpdateTrialData(hObject,eventData,handles)

parent = getappdata(handles.figure1,'parent');

trialData = getappdata(parent,'trialData');
currentDataSet = getappdata(parent,'CurrentDataSet');

expEventData = getappdata(parent,'expEventData');

if isempty(expEventData)
    error('No data set is loaded yet!')
end

evt = expEventData(currentDataSet); 

if ~isfield(evt,'events')  
    return
end
  
evttype = [expEventData(currentDataSet).events.type];
trialOnsets = double([expEventData(currentDataSet).events(evttype == 1).time]);
trialOnsetCodes = [expEventData(currentDataSet).events(evttype == 1).xdatcode];
trialOnsetCodeLabels = {expEventData(currentDataSet).events(evttype == 1).label};
trialOnsetEvtCode = {expEventData(currentDataSet).events(evttype == 1).code};

trialEnds = double([expEventData(currentDataSet).events(evttype == 2).time]);
trialEndCodes = [expEventData(currentDataSet).events(evttype == 2).xdatcode];
trialEndCodeLabels = {expEventData(currentDataSet).events(evttype == 2).label};
trialEndEvtCode = {expEventData(currentDataSet).events(evttype == 2).code};

if  length(trialOnsets)~=length(trialEnds)
    warning('Number of trialOnsets does not match the number of trial ends')
    return
elseif any(trialEnds-trialOnsets < 0)
    warning('Some trials have negative duration. Trial beginnings and ends are likely mismatched.')
end

if currentDataSet > length(trialData)
%     trialData(currentDataSet).codeincr = 0;
     if isempty(trialData)
         clear trialData;
     end
     trialData(currentDataSet) = makeTrialData([]);
end

if ~isempty(trialOnsets)    
    trialData(currentDataSet) = makeTrialData('startTime',trialOnsets,'stopTime', trialEnds,'startCode',trialOnsetCodes,...
        'stopCode',trialEndCodes,'code',(1:length(trialOnsets))+ trialData(currentDataSet).codeincr,...
        'startCodeLabel',trialOnsetCodeLabels,'stopCodeLabel',trialEndCodeLabels,'startEventCode',trialOnsetEvtCode,...
        'stopEventCode',trialEndEvtCode);

else
    trialData(currentDataSet) = makeTrialData;
end

% trialData(currentDataSet).trials = trialdat; 
% trialData(currentDataSet).codeincr = max([trialdat.code]);

setappdata(parent,'trialData',trialData);

tmfuns = getappdata(parent,'trialManagerFunctions');
tmhandle = getappdata(parent,'trialManager');

if ishandle(tmhandle)
    tmfuns.update;
end


% --- Executes on selection change in eventTypeMenu.
function eventTypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to eventTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns eventTypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eventTypeMenu


% --- Executes during object creation, after setting all properties.
function eventTypeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in markAllCodesCheck.
function markAllCodesCheck_Callback(hObject, eventdata, handles)
% hObject    handle to markAllCodesCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of markAllCodesCheck




% --------------------------------------------------------------------
function removeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to removeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

expEventData = getappdata(parent,'expEventData');
currentDataSet = getappdata(parent,'CurrentDataSet');

selected = get(handles.eventList,'value')-1;

selected(selected==0) = [];

expEventData(currentDataSet).events(selected) = [];

setappdata(parent,'expEventData',expEventData)

set(handles.eventList,'value',selected(1));
Update(hObject, eventdata, handles);


% --------------------------------------------------------------------
function eventListContextMenu_Callback(hObject, eventdata, handles)
% hObject    handle to eventListContextMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');
expEventData = getappdata(parent,'expEventData');
currentDataSet = getappdata(parent,'CurrentDataSet');
xdat = expEventData(currentDataSet).xdat;

xdatTs = [xdat.startT];
xdatcodenum = [xdat.id];
type = zeros(1,length(xdatTs));
xdatlabel = {xdat.code};


expEventData(currentDataSet) = makeEventData('time',xdatTs,'type',type,'xdatcode',xdatcodenum,'label',xdatlabel);
expEventData(currentDataSet).xdat = xdat;

setappdata(parent,'expEventData',expEventData);

UpdateTrialData(hObject,eventdata,handles)
Update(hObject,eventdata,handles)


% --- Executes on button press in duplicateCheckBox.
function duplicateCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to duplicateCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of duplicateCheckBox


