function varargout = GridManager(varargin)
% GRIDMANAGER M-file for GridManager.fig
%      GRIDMANAGER, by itself, creates a new GRIDMANAGER or raises the existing
%      singleton*.
%
%      H = GRIDMANAGER returns the handle to a new GRIDMANAGER or the handle to
%      the existing singleton*.
%
%      GRIDMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRIDMANAGER.M with the given input arguments.
%
%      GRIDMANAGER('Property','Value',...) creates a new GRIDMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GridManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GridManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GridManager

% Last Modified by GUIDE v2.5 16-Dec-2007 15:28:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GridManager_OpeningFcn, ...
                   'gui_OutputFcn',  @GridManager_OutputFcn, ...
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


% --- Executes just before GridManager is made visible.
function GridManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GridManager (see VARARGIN)

% Default Screen resolution, size and distance

defaultN = [15 15]; % Default number of divisions 
defaultdata = [0 1 0 1]; % Default range of grid (entire screen)
defaultLabel = 'Sampling_Grid';

parent = varargin{1};

setappdata(handles.figure1,'parent', parent );
setappdata(parent,'gridManager',hObject)
% gridmanagerfcns.updateBinData= @updateBinData;
gridmanagerfcns.updatePlot = @(varargin) updatePlot([],[],handles,varargin{:});
gridmanagerfcns.updateFields= @(varargin)updateFields([],[],handles,varargin{:});
setappdata(varargin{1},'gridManagerFunctions',gridmanagerfcns)

set(handles.figure1,'DeleteFcn',@Destructor)

if ~isappdata(parent,'binData')
    binData = makeBinData({defaultdata,defaultN},'type','grid','label',defaultLabel);
    currentBinGroup = 1;
    
else
%     currentGroup = getappdata(parent,'CurrentBinGroup');
%     if isempty(binData.groups(currentGroup).pos)
%          setappdata(parent,'CurrentBinGroup',currentGroup)
%     else
        binData = getappdata(parent,'binData');
    %     ngroups = length(binData.groups);
        types = {binData.groups.type};
        grid = find(strcmp(types,'grid'),1);
        if ~isempty(grid)
            currentBinGroup= grid(find(grid,1));
        else 
            currentBinGroup= length(binData.groups)+1 ;
            binData = makeBinData(binData,[] ,'type','grid');

%             binData.groups(currentBinGroup).pos = makeBinData([],'type','grid');
        end
        
%     end
end

setappdata(parent,'CurrentBinGroup',currentBinGroup);    
setappdata(parent,'binData',binData);

if ~ishandle(getappdata(parent,'binManager'))
    GazeReader('binManagerMenu_Callback',parent,[],guidata(parent));
end
    
    
       
%Sets the activecontrol flag in the parent object to 'gridManager'
MakeActive(hObject, eventdata, handles);

%Calls a function which assigns all fields to the value displayed in he
%gridManager window
% currentGroup = getappdata(parent,'CurrentBinGroup');
% if currentGroup > length(binData.groups) || isempty(binData.groups(currentGroup).pos)
%     binData.groups(currentGroup).inputData{1} = defaultdata;
%     binData.groups(currentGroup).inputData{2} = defaultN;    
%     binData.groups(currentGroup).label = defaultLabel;    
%    
% %     assignAllFields(hObject,eventdata,handles);
% end


if ~isempty(binData.groups) && ~isempty(binData.groups(currentBinGroup).pos)
    updateFields(hObject,eventdata,handles);
end
    updateBinData(hObject,eventdata,handles);
% updateFields(hObject,eventdata,handles);


% set(hObject,'WindowButtonDownFcn',  @FigureButtonDownFcn )
% set(handles.figure1,'WindowButtonUpFcn',  @FigureButtonUpFcn )

% Choose default command line output for GridManager
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GridManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GridManager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% gridmanagerfcns.updateBinData= @updateBinData;
% varargout{2} = gridmanagerfcns;

function xrange1_Callback(hObject, eventdata, handles)
% hObject    handle to xrange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MakeActive(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of xrange1 as text
%        str2double(get(hObject,'String')) returns contents of xrange1 as a double




updateBinData(hObject,eventdata,handles);
activefigure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function xrange1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xrange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xrange2_Callback(hObject, eventdata, handles)
% hObject    handle to xrange2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MakeActive(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of xrange2 as text
%        str2double(get(hObject,'String')) returns contents of xrange2 as a double



updateBinData(hObject,eventdata,handles);
activefigure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function xrange2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xrange2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yrange1_Callback(hObject, eventdata, handles)
% hObject    handle to yrange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MakeActive(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of yrange1 as text
%        str2double(get(hObject,'String')) returns contents of yrange1 as a double




updateBinData(hObject,eventdata,handles);
activefigure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function yrange1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yrange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yrange2_Callback(hObject, eventdata, handles)
% hObject    handle to yrange2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of yrange2 as text
%        str2double(get(hObject,'String')) returns contents of yrange2 as a double


updateBinData(hObject,eventdata,handles);
activefigure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function yrange2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yrange2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nx_Callback(hObject, eventdata, handles)
% hObject    handle to Nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of Nx as text
%        str2double(get(hObject,'String')) returns contents of Nx as a double




% setappdata(handles.figure1,'Nx',str2num(get(handles.Nx,'string')));

updateBinData(hObject, eventdata, handles)
activefigure(handles.figure1)


% --- Executes during object creation, after setting all properties.

function Nx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% setappdata(handles.figure1,'Ny',str2num(get(handles.Nx,'string')));


%---------------
% function grid = rect2grid(hObject, eventdata, handles,rect)
% 
% %Converts the input data format for rect type to that for grid
% 
% if iscell(rect)
%     grid = rect;
%     return
% end
% 
% minmax = cat(1,min(rect,[],1),max(rect,[],2));
% 
% grid{1} = minmax([1 4 5 8]);
% 
% nbins = str2num([get(handles.Nx,'string'),'  ',get(handles.Ny,'string')]);
% grid{2} = nbins;

%%%%%%%%%%%%%%%%%%%%%%

function Ny_Callback(hObject, eventdata, handles)
% hObject    handle to Ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ny as text
%        str2double(get(hObject,'String')) returns contents of Ny as a double

MakeActive(hObject, eventdata, handles)



% setappdata(handles.figure1,'Ny',str2num(get(handles.Nx,'string')));
updateBinData(hObject, eventdata, handles)
activefigure(handles.figure1)



% --- Executes during object creation, after setting all properties.
function Ny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% MakeActive(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%------------------------------------------------
  
function FigureButtonDownFcn(hObject,eventdata)

MakeActive(hObject, eventdata, handles)

setappdata(getappdata(hObject,'mainWin'),'activeWin','gridmanager');
MakeActive(hObject, eventdata, handles)

    


%-----------------------------------------------

% function updatePlot(hObject,eventdata,roiInstance,pivot,btnDownPt,axh)
function updateBinData(hObject,eventdata,handles)

parent = getappdata(handles.figure1,'parent');
binData = getappdata(parent,'binData');
currentGroup = getappdata(parent,'CurrentBinGroup');
screenData = getappdata(parent,'screenData')';
bmh = getappdata(parent,'binManager');
if isempty(bmh) || ~ishandle(bmh) % Open the bin manager window if not already
    GazeReader('binManagerMenu_Callback',parent,[],guidata(parent));
    MakeActive(hObject, eventdata, handles)
    activefigure(handles.figure1)
end

CurrentDataSet = getappdata(parent,'CurrentDataSet');
if isempty(CurrentDataSet) || CurrentDataSet == 0, CurrentDataSet = 1; end


if currentGroup == 0 || ~strcmp(binData.groups(currentGroup).type,'grid')
    return    
end

% data = binData.groups(currentGroup).inputData;
% data = binData.groups(currentGroup).pos;
% data = rect2grid(data);

data{2}(1) = str2double(get(handles.Nx,'String'));
data{2}(2) = str2double(get(handles.Ny,'String'));
data{1}(1) = str2double(get(handles.xrange1,'String'))./screenData(CurrentDataSet).res(1);
data{1}(2) = str2double(get(handles.xrange2,'String'))./screenData(CurrentDataSet).res(1);
data{1}(3) = str2double(get(handles.yrange1,'String'))./screenData(CurrentDataSet).res(2);
data{1}(4) = str2double(get(handles.yrange2,'String'))./screenData(CurrentDataSet).res(2);


% Create a new bin group with updated bin position data

label = binData.groups(currentGroup).label;
code = binData.groups(currentGroup).code;

newBinData = makeBinData( data ,'trials',binData.groups(currentGroup).activeTrials,'precedence',binData.precedence,'label',label,'type','grid','code',code);
binData.groups(currentGroup)= newBinData.groups;
    
setappdata(parent,'binData',binData);

updateFields(hObject,eventdata,handles)
updatePlot(hObject,eventdata,handles);

%-----------------------------------------------

% function updatePlot(hObject,eventdata,roiInstance,pivot,btnDownPt,axh)
function updatePlot(hObject,eventdata,handles, varargin)


% plottype = 'g'; %arguments to plot for drawing bin lines

% Simply calls the plotting function handle associated with this bin group

if ~ishandle(handles.figure1)
    parent = varargin{1};
else
    parent = getappdata(handles.figure1,'parent');
end

if nargin > 3 && ishandle(varargin{1})
    varargin(1) = [];
end

% phandles = guidata(parent);
% % axes(phandles.axes2)
% binData = getappdata(parent,'binData');
binManagerFunctions = getappdata(parent,'binManagerFunctions');

setappdata(parent,'CurrentBin',0);

if ~isempty(binManagerFunctions)
    binManagerFunctions.draw();
end
% currim = getappdata(parent,'CurrentImage');
% trialData = getappdata(parent,'trialData');
% screenData = getappdata(parent,'screenData');


% currentGroup = getappdata(parent,'CurrentBinGroup');
% 
% 
% plothandles = getappdata(handles.figure1,'plotHandles');
% 
% delete(plothandles(ishandle(plothandles(:))));
% 
% % axes(phandles.axes2)
% 
% % axis(phandles.axes2,axis(phandles.axes1));
% ax1 = get(phandles.axes1);
% set(phandles.axes2,'units',ax1.Units,'position',ax1.Position,'Ydir',ax1.YDir)
% 
% hold(phandles.axes2,'off')
% plothandles = binData.groups(currentGroup).plot(binData.groups(currentGroup),0,plottype,'parent',phandles.axes2, varargin{:});
% % hold off
% axis(phandles.axes2,'off')
% % axes(phandles.axes1)
% 
% setappdata(handles.figure1,'plotHandles',plothandles)
% 

%---------------

function Destructor(hObject,eventdata,handles)

ph = getappdata(hObject,'plotHandles');
delete(ph(ishandle(ph)))

%-----------------------------------------------

% function updatePlot(hObject,eventdata,roiInstance,pivot,btnDownPt,axh)
function updateFields(hObject,eventdata,handles, varargin)
    
%Sets all fields in the gridManager to the current settings contained in binData


if ~ishandle(handles.figure1)
    parent = varargin{1};
else
    parent = getappdata(handles.figure1,'parent');
end


binData =  getappdata(parent,'binData');

CurrentBinGroup = getappdata(parent,'CurrentBinGroup');

CurrentDataSet = getappdata(parent,'CurrentDataSet');
if isempty(CurrentDataSet) || CurrentDataSet == 0, CurrentDataSet = 1; end


data = binData.groups(CurrentBinGroup).pos;
if ~iscell(data)
    data = rect2grid(data);
end
label= binData.groups(CurrentBinGroup).label;
screenData = getappdata(parent,'screenData');
res = screenData(CurrentDataSet).res;

set(handles.xrange1,'String',num2str(data{1}(1)*res(1) ) );
set(handles.xrange2,'String',num2str(data{1}(2)*res(1) ) );
set(handles.yrange1,'String',num2str(data{1}(3)*res(2) ) );
set(handles.yrange2,'String',num2str(data{1}(4)*res(2) ) );
set(handles.Nx,'String',num2str(data{2}(1)));
set(handles.Ny,'String',num2str(data{2}(2)));
set(handles.label,'String',label);



%------------------------------------------------
  
function FigureButtonUpFcn(hObject,eventdata)

% setappdata(hObject,'buttonState',0)
% setappdata(hObject,'buttonNumber',0)
% 
% MakeActive(hObject, eventdata, handles)

% --- Executes on button press in equalWidth.
function equalWidth_Callback(hObject, eventdata, handles)
% hObject    handle to equalWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of equalWidth

MakeActive(hObject, eventdata, handles)
activefigure(handles.figure1)

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
activefigure(handles.figure1)

% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
activefigure(handles.figure1)

%%%%%%%%%%%%%%%%%%%%%%
function MakeActive(hObject, eventdata, handles)
% Makes the current control active in the main window

parent = getappdata(handles.figure1,'parent');

setappdata(parent,'activeControl','gridManager');



function label_Callback(hObject, eventdata, handles)
% hObject    handle to label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label as text
%        str2double(get(hObject,'String')) returns contents of label as a double


MakeActive(hObject, eventdata, handles)


parent = getappdata(handles.figure1,'parent');
binData = getappdata(parent,'binData');
currentGroup = getappdata(parent,'CurrentBinGroup');
if currentGroup == 0
    return
end

binData.groups(currentGroup).label = get(handles.label,'String');
setappdata(parent,'binData',binData);
activefigure(handles.figure1)

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

%-------------------------------------

% function assignAllFields(hObject, eventdata, handles)
% 
% %Updates data in binData based on the current entry in all of the fields
% %of the gridManagerWindow
% 
% 
% parent = getappdata(handles.figure1,'parent');
% binData = getappdata(parent,'binData');
% currentGroup = getappdata(parent,'CurrentBinGroup');
% screenData = getappdata(parent,'screenData');
% res = screenData.res;
% 
% binData.groups(currentGroup).label = get(handles.label,'String');
% binData.groups(currentGroup).inputData{1} = str2num(get(handles.xrange1,'String'))./res(1);
% binData.groups(currentGroup).inputData{1}(2) = str2num(get(handles.xrange2,'String'))./res(1);
% binData.groups(currentGroup).inputData{1}(3) = str2num(get(handles.yrange1,'String'))./res(2);
% binData.groups(currentGroup).inputData{1}(4) = str2num(get(handles.yrange2,'String'))./res(2);
% binData.groups(currentGroup).inputData{2} = str2num(get(handles.Nx,'String'));
% binData.groups(currentGroup).inputData{2}(2) = str2num(get(handles.Ny,'String'));
% 
% setappdata(parent,'binData',binData);
% 
% updateBinData(hObject, eventdata, handles)
% 


