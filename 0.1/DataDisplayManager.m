function varargout = DataDisplayManager(varargin)
% DATADISPLAYMANAGER M-file for DataDisplayManager.fig
%      DATADISPLAYMANAGER, by itself, creates a new DATADISPLAYMANAGER or raises the existing
%      singleton*.
%
%      H = DATADISPLAYMANAGER returns the handle to a new DATADISPLAYMANAGER or the handle to
%      the existing singleton*.
%
%      DATADISPLAYMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATADISPLAYMANAGER.M with the given input arguments.
%
%      DATADISPLAYMANAGER('Property','Value',...) creates a new DATADISPLAYMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DataDisplayManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DataDisplayManager_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help DataDisplayManager

% Last Modified by GUIDE v2.5 12-Jan-2008 01:55:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataDisplayManager_OpeningFcn, ...
                   'gui_OutputFcn',  @DataDisplayManager_OutputFcn, ...
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


% --- Executes just before DataDisplayManager is made visible.
function DataDisplayManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataDisplayManager (see VARARGIN)

parent = varargin{1};

dmfun.update = @(varagin)Update([],[],handles);
dmfun.drawRaw= @(varagin)DrawRaw([],[],handles,parent,varargin{:});
dmfun.drawFix= @(varagin)DrawFix([],[],handles,parent,varargin{:});


setappdata(handles.figure1,'parent',parent);
setappdata(parent,'DataDisplayFunctions',dmfun);
setappdata(parent,'DataDisplayManager',handles.figure1);


% Choose default command line output for DataDisplayManager
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

Update(hObject, eventdata, handles)
% UIWAIT makes DataDisplayManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DataDisplayManager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in fixationList.
function fixationList_Callback(hObject, eventdata, handles)
% hObject    handle to fixationList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fixationList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fixationList

parent = getappdata(handles.figure1,'parent');

selected = get(handles.fixationList,'value');
CurrentDataSet = getappdata(parent,'CurrentDataSet');
trialData = getappdata(parent,'trialData');
currentTrial = getappdata(parent,'CurrentTrial');


if selected == 1 || currentTrial == 0
    
    setappdata(parent,'CurrentFixation',selected-1)

else
    
    setappdata(parent,'CurrentFixation',trialData(CurrentDataSet).trials(currentTrial).fixations(selected-1))
end

Update(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function fixationList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixationList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fixEdit_Callback(hObject, eventdata, handles)
% hObject    handle to fixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixEdit as text
%        str2double(get(hObject,'String')) returns contents of fixEdit as a double
parent = getappdata(handles.figure1,'parent');

CurrentDataSet = getappdata(parent,'CurrentDataSet');
trialData = getappdata(parent,'trialData');


newfx = str2double(get(handles.fixEdit,'string'));

currtd = trialData(CurrentDataSet);
ntrials = length(currtd.trials);



if ~isempty(newfx) 
    
      tr = find(cellfun(@ismember,num2cell(newfx*ones(1,ntrials)),{trialData(CurrentDataSet).trials.fixations})); 
    setappdata(parent,'CurrentFixation',newfx );
    setappdata(parent,'CurrentTrial',tr);
end

Update(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function fixEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trialEdit_Callback(hObject, eventdata, handles)
% hObject    handle to trialEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trialEdit as text
%        str2double(get(hObject,'String')) returns contents of trialEdit as a double

parent = getappdata(handles.figure1,'parent');

CurrentDataSet = getappdata(parent,'CurrentDataSet');
trialData = getappdata(parent,'trialData');


newtr = str2double(get(handles.trialEdit,'string'));

if ~isempty(newtr) && newtr <= length(trialData(CurrentDataSet).trials)

    setappdata(parent,'CurrentTrial',newtr);
end
Update(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function trialEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in jumpFixUp.
function jumpFixUp_Callback(hObject, eventdata, handles)
% hObject    handle to jumpFixUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');

currentTrial = getappdata(parent,'CurrentTrial');
currentFixation = getappdata(parent,'CurrentFixation');
CurrentDataSet = getappdata(parent,'CurrentDataSet');

trialData = getappdata(parent,'trialData');


if currentTrial==0
    currentTrial = currentTrial+1;
    setappdata(parent,'CurrentTrial',currentTrial);
    setappdata(parent,'CurrentFixation',trialData(CurrentDataSet).trials(currentTrial).fixations(1));

else
    tr = trialData(CurrentDataSet).trials(currentTrial);

    if tr.fixations(end) == currentFixation

        if currentTrial == length(trialData(CurrentDataSet).trials)
            return
        end
        currentTrial = currentTrial+1;
        setappdata(parent,'CurrentTrial',currentTrial);
        setappdata(parent,'CurrentFixation',trialData(CurrentDataSet).trials(currentTrial).fixations(1));

    else
        setappdata(parent,'CurrentFixation',currentFixation+1);
    end
end
    

Update(hObject, eventdata, handles)

% --- Executes on button press in jumpTrialUp.
function jumpTrialUp_Callback(hObject, eventdata, handles)
% hObject    handle to jumpTrialUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');

currentTrial = getappdata(parent,'CurrentTrial');
CurrentDataSet = getappdata(parent,'CurrentDataSet');

trialData = getappdata(parent,'trialData');

if currentTrial < length(trialData(CurrentDataSet).trials)
    setappdata(parent,'CurrentTrial',currentTrial+1);
end
    
setappdata(parent,'CurrentFixation',0);

Update(hObject, eventdata, handles)


% --- Executes on button press in jumpTrialDown.
function jumpTrialDown_Callback(hObject, eventdata, handles)
% hObject    handle to jumpTrialDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

currentTrial = getappdata(parent,'CurrentTrial');

if currentTrial > 0
    setappdata(parent,'CurrentTrial',currentTrial-1);
end
    
setappdata(parent,'CurrentFixation',0);

Update(hObject, eventdata, handles)

% --- Executes on button press in jumpFixDown.
function jumpFixDown_Callback(hObject, eventdata, handles)
% hObject    handle to jumpFixDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');

currentTrial = getappdata(parent,'CurrentTrial');
currentFixation = getappdata(parent,'CurrentFixation');
CurrentDataSet = getappdata(parent,'CurrentDataSet');

trialData = getappdata(parent,'trialData');


if currentTrial > 0 && currentFixation == trialData(CurrentDataSet).trials(currentTrial).fixations(1)
    setappdata(parent,'CurrentTrial',currentTrial-1);
    
    if currentTrial > 1
        tr = trialData(CurrentDataSet).trials(currentTrial-1);
        if ~isempty(tr.fixations)
            setappdata(parent,'CurrentFixation',tr.fixations(end));
        end
    end
elseif currentFixation  > 0 
    
    setappdata(parent,'CurrentFixation',currentFixation-1);

end

Update(hObject, eventdata, handles)

%------------------------------------------
function Update(hObject, eventdata, handles)


parent = getappdata(handles.figure1,'parent');


currentTrial = getappdata(parent,'CurrentTrial');

currentFixation = getappdata(parent,'CurrentFixation');
if isempty(currentFixation)
    currentFixation = 0;
end

CurrentDataSet = getappdata(parent,'CurrentDataSet');
trialData = getappdata(parent,'trialData');
fixdata = getappdata(parent,'fixationData');
pld = getappdata(parent,'DataPlotHandles');

delete(pld(ishandle(pld)));


if currentTrial~= 0 && isempty(trialData(CurrentDataSet).trials(currentTrial).startTime)
    return
end
if get(handles.rawDisplayButton,'value')

    DrawRaw(hObject,eventdata,handles,parent,'y')

    
end        
        
if get(handles.fixationDisplay,'value')

    DrawFix(hObject,eventdata,handles,parent,'w')

    
end        
set(handles.trialEdit,'string',num2str(currentTrial))
set(handles.fixEdit,'string',num2str(currentFixation))

if get(handles.fixationDisplay,'value') == 0
    fxs = [];
elseif  currentTrial == 0  
    
    fxs = 1:length(fixdata(CurrentDataSet).fix);    
elseif currentTrial > 0
    fxs = trialData(CurrentDataSet).trials(currentTrial).fixations;
end


liststr = {'#      Start Time      Duration'};

if ~isempty(fxs)
    fixnum = cellfun(@num2str,mat2cell(fxs,1,ones(1,length(fxs))),'uniformoutput',0);
    cstrlen = cellfun(@length,fixnum);
    spaces1 = cellfun(@(n) repmat(' ',1,n),mat2cell(15-2*cstrlen,1,ones(1,length(cstrlen))),'UniformOutput',false);

    fixstt = cellfun(@num2str,{fixdata(CurrentDataSet).fix(fxs).startT},'uniformoutput',0);
    cstrlen = cellfun(@length,fixstt);
    spaces2 = cellfun(@(n) repmat(' ',1,n),mat2cell(15-2*cstrlen,1,ones(1,length(cstrlen))),'UniformOutput',false);

    fixdur = cellfun(@num2str,{fixdata(CurrentDataSet).fix(fxs).dur},'uniformoutput',0);
    cstrlen = cellfun(@length,fixdur);
    spaces3 = cellfun(@(n) repmat(' ',1,n),mat2cell(15-2*cstrlen,1,ones(1,length(cstrlen))),'UniformOutput',false);

    liststr = cat(2,liststr ,strcat(fixnum,spaces1,fixstt,spaces2,fixdur,spaces3));

    set(handles.fixationList,'string',liststr)
end
if currentTrial == 0 || currentFixation == 0
    set(handles.fixationList,'value',currentFixation+1)
else
    set(handles.fixationList,'value',find(fxs == currentFixation)+1)
end    

rmfun = getappdata(parent,'RegManagerFunctions');
bmfun = getappdata(parent,'binManagerFunctions');
trfun = getappdata(parent,'trialManagerFunctions');

if ~isempty(bmfun) && get(handles.fixationDisplay,'value') 
    bmfun.update(parent,false);
end

if ~isempty(rmfun)
    rmfun.update();
end
if ~isempty(trfun)
    trfun.update();
end



%------------------
function DrawRaw(hObject, eventdata, handles,parent,varargin)

if ~isempty(parent) || nargin < 5 
    parent = getappdata(handles.figure1,'parent');
end


trialData = getappdata(parent,'trialData');

currentTrial = getappdata(parent,'CurrentTrial');

currentFixation = getappdata(parent,'CurrentFixation');
CurrentDataSet = getappdata(parent,'CurrentDataSet');



scrdat = getappdata(parent,'screenData');

phandles = guidata(parent);

rawdata = getappdata(parent,'rawGazeData');
    
 disp = 1:length(rawdata(CurrentDataSet).horz);
 time = rawdata(CurrentDataSet).time;
if currentTrial ~= 0
  
    disp =time >=  trialData(CurrentDataSet).trials(currentTrial).startTime & time  < trialData(CurrentDataSet).trials(currentTrial).stopTime;

%     disp = disp - rawdata(CurrentDataSet).time(1)+1;
end

pld = getappdata(parent,'DataPlotHandles');

hold(phandles.axes2,'on')

pld = cat(1,pld,plot(rawdata(CurrentDataSet).horz(disp)./scrdat(CurrentDataSet).res(1),rawdata(CurrentDataSet).vert(disp)./scrdat(CurrentDataSet).res(2),varargin{:},'parent',phandles.axes2));


if  currentFixation ~=0
    fixdata= getappdata(parent,'fixationData');
    
    disp =time >=  fixdata(CurrentDataSet).fix(currentFixation).startT  & time  < fixdata(CurrentDataSet).fix(currentFixation).startT+fixdata(CurrentDataSet).fix(currentFixation).dur;

%     disp = fixdata(CurrentDataSet).fix(currentFixation).startT + (0:fixdata(CurrentDataSet).fix(currentFixation).dur-1);
%     disp = disp - rawdata(CurrentDataSet).time(1)+1;
    pld = cat(1,pld,plot(rawdata(CurrentDataSet).horz(disp)./scrdat(CurrentDataSet).res(1),rawdata(CurrentDataSet).vert(disp)./scrdat(CurrentDataSet).res(2),'r','parent',phandles.axes2));
end
setappdata(parent,'DataPlotHandles',pld);


%------------------
function DrawFix(hObject, eventdata, handles,parent,varargin)

CurrentDataSet= getappdata(parent,'CurrentDataSet');


fixationData = getappdata(parent,'fixationData');

% durscale =1./10000; %  Scaling of Duration circle radius
durscale = .01./std(double([fixationData(CurrentDataSet).fix.dur]));

posField = 'meanPos';
if ~isempty(parent) || nargin < 5 
    parent = getappdata(handles.figure1,'parent');
end

trialData = getappdata(parent,'trialData');

currentTrial = getappdata(parent,'CurrentTrial');

currentFixation = getappdata(parent,'CurrentFixation');

pld = getappdata(parent,'DataPlotHandles');

scrdat = getappdata(parent,'screenData');

phandles = guidata(parent);

if  currentTrial == 0 && currentFixation == 0        
    fxs = 1:length(fixationData(CurrentDataSet).fix);
elseif  currentTrial == 0       
    fxs = currentFixation ;
else
    fxs = trialData(CurrentDataSet).trials(currentTrial).fixations;
end

fixposraw = cat(1,fixationData(CurrentDataSet).fix(fxs).(posField));
if isempty(fixposraw), return, end
fxpos = (fixposraw*diag(scrdat(CurrentDataSet).res.^-1))';
fxdur = double(cat(2,fixationData(CurrentDataSet).fix(fxs).dur))*durscale;

hold(phandles.axes2,'on')

pld = cat(1,pld,plot(fxpos(1,:),fxpos(2,:),'parent',phandles.axes2));

th = (0:.1:2*pi)';

circx = cos(th)*fxdur + repmat(fxpos(1,:),length(th),1);
circy = sin(th)*fxdur + repmat(fxpos(2,:),length(th),1);

pld = cat(1,pld,plot(circx,circy,varargin{:},'parent',phandles.axes2));

if ~isempty(currentFixation) && currentFixation ~= 0 && ismember(currentFixation,fxs)
    pld = cat(1,pld,plot(circx(:,ismember(fxs,currentFixation)),circy(:,ismember(fxs,currentFixation)),'r','parent',phandles.axes2));
    pld = cat(1,pld,plot(fxpos(1,ismember(fxs,currentFixation)),fxpos(2,ismember(fxs,currentFixation)),'r','parent',phandles.axes2));
end

setappdata(parent,'DataPlotHandles',pld);


% --- Executes on button press in fixationDisplay.
function fixationDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to fixationDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixationDisplay


Update(hObject, eventdata, handles)


% --- Executes on button press in rawDisplayButton.
function rawDisplayButton_Callback(hObject, eventdata, handles)
% hObject    handle to rawDisplayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rawDisplayButton


Update(hObject, eventdata, handles)


