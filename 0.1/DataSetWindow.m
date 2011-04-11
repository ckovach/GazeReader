function varargout = DataSetWindow(varargin)
% DATASETWINDOW M-file for DataSetWindow.fig
%      DATASETWINDOW, by itself, creates a new DATASETWINDOW or raises the existing
%      singleton*.
%
%      H = DATASETWINDOW returns the handle to a new DATASETWINDOW or the handle to
%      the existing singleton*.
%
%      DATASETWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATASETWINDOW.M with the given input arguments.
%
%      DATASETWINDOW('Property','Value',...) creates a new DATASETWINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DataSetWindow_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DataSetWindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DataSetWindow

% Last Modified by GUIDE v2.5 21-Apr-2008 11:06:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataSetWindow_OpeningFcn, ...
                   'gui_OutputFcn',  @DataSetWindow_OutputFcn, ...
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


% --- Executes just before DataSetWindow is made visible.
function DataSetWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataSetWindow (see VARARGIN)

% Choose default command line output for DataSetWindow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

parent = varargin{1};
setappdata(handles.figure1,'parent',parent);


dswfun.pool = @(datasets)Pool_Callback(hObject,[],handles,datasets);
dswfun.delete = @(varargin)removeDataSet_Callback(hObject,[],handles,varargin{:});

setappdata(parent,'DataSetWindowFunctions',dswfun);

UpdateFields(hObject,eventdata,handles);




% UIWAIT makes DataSetWindow wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DataSetWindow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function label_Callback(hObject, eventdata, handles)
% hObject    handle to label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label as text
%        str2double(get(hObject,'String')) returns contents of label as a double

parent = getappdata(handles.figure1,'parent');

datainfo = getappdata(parent,'eyetrackerHeaderData');
currentDataSet = getappdata(parent,'CurrentDataSet');

if length(datainfo) == 1
    currentDataSet = 1;
    setappdata(parent,'eyetrackerHeaderData',currentDataSet);
elseif currentDataSet == 0
    return
end

datainfo(currentDataSet).label = get(handles.label,'string');

setappdata(parent,'eyetrackerHeaderData',datainfo);

UpdateFields(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dataList.
function dataList_Callback(hObject, eventdata, handles)
% hObject    handle to dataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dataList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataList

parent = getappdata(handles.figure1,'parent');


value = get(handles.dataList,'value');
if length(value) == 1
    setappdata(parent,'CurrentDataSet',value-1);
    UpdateFields(hObject, eventdata, handles)
end




% --- Executes during object creation, after setting all properties.
function dataList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%

function UpdateFields(hObject,eventData,handles)

%Updates the data set window 

parent = getappdata(handles.figure1,'parent');
currentDataSet = getappdata(parent,'CurrentDataSet');

if isempty(currentDataSet) 
    set(handles.label,'String',[]);
    set(handles.filename,'String','None Selected');    
    set(handles.dataList,'value',1);
    return
end    

datainfo = getappdata(parent,'eyetrackerHeaderData');

if currentDataSet > length(datainfo) 
    currentDataSet = 0;
end
if ~isempty(datainfo)
    filenames = {datainfo.filename};
%     [path,filename] = fileparts(filenames{currentDataSet });

    set(handles.filename,'String',filenames);

    setnums = cellfun(@num2str,mat2cell(1:length(datainfo),1,ones(1,length(datainfo))),'uniformOutput',0);

    labels = {datainfo.label};
    if currentDataSet~=0
        set(handles.label,'String',labels{currentDataSet});
    else
        set(handles.label,'String','none selected');
    end
else
    currentDataSet = 0;
    setappdata(parent,'CurrentDataSet' ,currentDataSet );
    setnums={};
    labels= {};
    filenames = {};
end
liststr = [{'#           Label               File'},strcat(setnums,{'       '},labels,{'       '},filenames)];
set(handles.dataList,'String',liststr);
set(handles.dataList,'value',currentDataSet+1);

if currentDataSet > 0
    set(handles.filename,'string',datainfo(currentDataSet).filename);
    set(handles.label,'string',datainfo(currentDataSet).label);
else
    set(handles.filename,'string','none selected');
    set(handles.label,'string','none selected');
end


% --------------------------------------------------------------------
function removeDataSet_Callback(hObject, eventdata, handles,varargin)
% hObject    handle to removeDataSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

if isempty(varargin)
    DataSets = get(handles.dataList,'value')-1;
else
    DataSets = varargin{1};
end
DataSets = DataSets(DataSets~=0);


datainfo = getappdata(parent,'eyetrackerHeaderData');
% data= getappdata(parent,'eyetrackerData');

fixationData = getappdata(parent,'fixationData');
rawData = getappdata(parent,'rawGazeData');
expEventData = getappdata(parent,'expEventData');
regData = getappdata(parent,'regData');
modelData = getappdata(parent,'modelData');
screenData = getappdata(parent,'screenData');
try
    fixationData(DataSets) = [];
catch
    warning('Size of fixation Data  data structure doesn''t match fixation data')
end    
try
    rawData(DataSets) = [];
catch
    warning('Size of raw  data structure doesn''t match fixation data')
   
end
try
    screenData(DataSets) = [];
catch
    warning('Size of raw  data structure doesn''t match fixation data')
   
end
try
   expEventData(DataSets) = [];
catch
   
    warning('Size of raw  data structure doesn''t match fixation data')
end
trialData = getappdata(parent,'trialData');
if ~isempty(trialData) 
    trialData(DataSets(DataSets <= length(trialData))) = [];
end
setappdata(parent,'trialData',trialData)

if ~isempty(regData) 
    regData(DataSets(DataSets <= length(regData))) = [];
end
setappdata(parent,'regData',regData)

if ~isempty(modelData) 
    modelData(DataSets(DataSets <= length(modelData))) = [];
end
setappdata(parent,'modelData',modelData)

% %Sets conditional model parameter 
% ConditionalModel = getappdata(parent ,'ConditionalModel');
% ConditionalModel(DataSets(DataSets <= length(modelData))) = [];
% TSampIntvl = getappdata(parent ,'TSampIntvl'); % sampling interval for binnin in time
% TSampIntvl(DataSets(DataSets <= length(modelData))) = [];
% etappdata(parent ,'ConditionalModel',ConditionalModel);
% setappdata(parent ,'TSampIntvl',TSampIntvl); % sampling interval for binnin in time





datainfo(DataSets) = [];
% data(DataSets) = [];

setappdata(parent,'fixationData',fixationData);
% setappdata(parent,'eyetrackerData',data);
setappdata(parent,'rawGazeData',rawData);
setappdata(parent,'expEventData',expEventData);

setappdata(parent,'eyetrackerHeaderData',datainfo);
setappdata(parent,'CurrentDataSet',DataSets(1)-1);
setappdata(parent,'CurrentTrial',0);
setappdata(parent,'CurrentFixation',0);
% setappdata(parent,'eyetrackerData',data);

UpdateFields(hObject,eventdata,handles);

    


% --------------------------------------------------------------------
function ListMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ListMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function Pool_Callback(hObject, eventdata, handles,varargin)
% hObject    handle to Pool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%concatenates data sets.

parent = getappdata(handles.figure1,'parent');
if isempty(varargin)
    DataSets = get(handles.dataList,'value')-1;
    DataSets = DataSets(DataSets~=0);
else
    DataSets = varargin{1};
end

regData = getappdata(parent,'regData');
% binData = getappdata(parent,'binData')
eventData = getappdata(parent,'expEventData');
rawGazeData = getappdata(parent,'rawGazeData');
trialData = getappdata(parent,'trialData');
fixationData = getappdata(parent,'fixationData');


[fixationData(end+1), lastendts, lastfxns,xdmaps] = concatFD(fixationData(DataSets));

[eventData(end+1),evtcodemaps] = concatED(lastendts, xdmaps,fixationData(end).xdatCodes,eventData(DataSets));

if   ~isempty(trialData) &&  min(DataSets) <=length(trialData)
    [trialData(end+1),lasttrs] = concatTD(lastendts,lastfxns,zeros(size(DataSets)),xdmaps, fixationData(end).xdatCodes,evtcodemaps,trialData(DataSets));
end
if ~isempty(regData) && min(DataSets) <=length(regData)
    regData(end+1) = concatRD(lastfxns,lasttrs,regData(DataSets));
end
if ~isempty(rawGazeData) && min(DataSets) <=length(rawGazeData)
    rawGazeData(end+1) = concatRWD(lastendts, rawGazeData(DataSets));
end

%Image Data !
 

datainfo = getappdata(parent,'eyetrackerHeaderData');
datainfo(end+1).label = ['pooled ',sprintf('%i ',DataSets)];
currentDataSet = length(datainfo);
setappdata(parent,'CurrentDataSet',currentDataSet);
setappdata(parent,'eyetrackerHeaderData',datainfo);
setappdata(parent,'fixationData',fixationData);
setappdata(parent,'regData',regData);
setappdata(parent,'trialData',trialData);
setappdata(parent,'rawGazeData',rawGazeData);
setappdata(parent,'expEventData',eventData);

%Sets conditional model parameter to that of the first data set
% ConditionalModel = getappdata(parent ,'ConditionalModel');
% TSampIntvl = getappdata(parent ,'TSampIntvl'); % sampling interval for binnin in time
% ConditionalModel(end+1) = ConditionalModel(DataSets(1));
% TSampIntvl(end+1) = TSampIntv(DataSets(1));
% setappdata(parent ,'ConditionalModel',ConditionalModel);
% setappdata(parent ,'TSampIntvl',TSampIntvl); % sampling interval for binnin in time
% 
% 
% 
