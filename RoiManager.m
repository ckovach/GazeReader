function varargout = RoiManager(varargin)
% ROIMANAGER M-file for RoiManager.fig
%      ROIMANAGER, by itself, creates a new ROIMANAGER or raises the existing
%      singleton*.
%
%      H = ROIMANAGER returns the handle to a new ROIMANAGER or the handle to
%      the existing singleton*.
%
%      ROIMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROIMANAGER.M with the given input
%      arguments.
%
%      ROIMANAGER('Property','Value',...) creates a new ROIMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RoiManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RoiManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RoiManager

% Last Modified by GUIDE v2.5 02-Jan-2008 17:34:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RoiManager_OpeningFcn, ...
                   'gui_OutputFcn',  @RoiManager_OutputFcn, ...
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


% --- Executes just before RoiManager is made visible.
function RoiManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RoiManager (see VARARGIN)

% Choose default command line output for RoiManager
handles.output = hObject;

parent =varargin{1};
setappdata(handles.figure1,'parent',parent );
MakeActive(hObject, eventdata, handles)

setappdata(parent,'CurrentRoi',0);
% currentRoiGroup = getappdata(parent,'CurrentRoiGroup');

roiManagerFunctions.draw = @(varargin)DrawPatch([],[],handles,varargin{:});
roiManagerFunctions.clear = @(varargin)ClearPatch([],[],handles,varargin{:});
roiManagerFunctions.update = @(varargin)UpdateRois([],[],handles,varargin{:});
% roiManagerFunctions.update = @SetFields;

setappdata(parent,'roiManagerFunctions',roiManagerFunctions);

currentRoiGroup = getappdata(parent,'CurrentRoiGroup');
roiData = getappdata(parent,'roiData');

if isempty(currentRoiGroup)
    setappdata(parent,'CurrentRoiGroup',0) ;
end

types = {'rect','ellipse','poly','grid'};
setappdata(handles.figure1,'types',types);

if isempty(roiData)
    
    roiData = makeBinData([],'label','Roi Group 1','type',types{ get( handles.type,'value' ) },'code',1 );
    setappdata(parent,'roiData',roiData)
    currentRoiGroup = 1;
%     roiData.codeincr = 0;
    setappdata(parent,'CurrentRoiGroup',currentRoiGroup)
    
end

UpdateRois(hObject, eventdata, handles)
% if currentRoiGroup ~= 0
%     SetFields(hObject, eventdata, handles);
% end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RoiManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RoiManager_OutputFcn(hObject, eventdata, handles) 
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
MakeActive(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');
currRoiGr = getappdata(parent,'CurrentRoiGroup');

if currRoiGr ~= 0 
    roiData = getappdata(parent,'roiData');
    roiData.groups(currRoiGr).label = get(handles.label,'String');
    setappdata(parent,'roiData',roiData)
end


SetFields(hObject, eventdata, handles)
UpdateRois(hObject, eventdata, handles)


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


% --- Executes on selection change in type.
function type_Callback(hObject, eventdata, handles)
% hObject    handle to type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type
MakeActive(hObject, eventdata, handles)

typeLabels = getappdata(handles.figure1,'types');

parent = getappdata(handles.figure1,'parent');
currRoiGr = getappdata(parent,'CurrentRoiGroup');

roiData = getappdata(parent,'roiData');
val = get(handles.type,'value');

if currRoiGr ~= 0 
%      strs = get(handles.type,'string');

    if ~isempty(roiData.groups(currRoiGr).pos)
        set(handles.type,'value',find(ismember(typeLabels,roiData.groups(currRoiGr).type)));
    else
        roiData.groups(currRoiGr).type = typeLabels{val};
        setappdata(parent,'roiData',roiData)
    end        
end

SetFields(hObject, eventdata, handles)
UpdateRois(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in roiGroupList.
function roiGroupList_Callback(hObject, eventdata, handles)
% hObject    handle to roiGroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns roiGroupList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiGroupList
MakeActive(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');

roiData = getappdata(parent,'roiData');

roival = get(handles.roiGroupList,'Value');

if length(roival) > 1
    ClearPatch(hObject, eventdata, handles);
    DrawPatch(hObject, eventdata, handles);

    return
end

setappdata(parent,'CurrentRoi',0);
currentBinGroup = roival-1;
setappdata(parent,'CurrentRoiGroup',currentBinGroup )    
if roival == 1 && strcmp(get(handles.figure1,'selectiontype'),'open')
    
    currRoiGr = length(roiData.groups)+1;
    
    setappdata(parent,'CurrentRoiGroup',currRoiGr);
%     setappdata(parent,'CurrentRoi',0);
    
    types = getappdata(handles.figure1,'types');
    type = types{get(handles.type,'Value')};
    roiData = makeBinData(roiData,[],'type',type,'label',sprintf('Roi Group %i',currRoiGr));
    
%     roiData.groups(currRoiGr) = newroi.groups;
    setappdata(parent,'roiData',roiData)
    
end

bgstr  = {roiData.groups.label};
set(handles.roiGroupList,'Value',currentBinGroup+1);
set(handles.roiGroupList,'String',cat(2,{'New ROI Group'},bgstr));

SetFields(hObject, eventdata, handles)
UpdateRois(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function roiGroupList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiGroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in roiList.
function roiList_Callback(hObject, eventdata, handles)
% hObject    handle to roiList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns roiList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiList


parent = getappdata(handles.figure1,'parent');
roiData = getappdata(parent,'roiData');

currentRoiGroup = getappdata(parent,'CurrentRoiGroup');
roinum = get(handles.roiList,'Value')-1;
types = getappdata(handles.figure1,'types');
type = types{get(handles.type,'value')};
label = get(handles.label,'string');

if length(roinum) > 1
    ClearPatch(hObject, eventdata, handles);
    DrawPatch(hObject, eventdata, handles);
    return
end

if currentRoiGroup ~=0 && ~isempty(roinum) && roinum  <= size(roiData.groups(currentRoiGroup).pos,1)
    setappdata(parent,'CurrentRoi',roinum)
elseif strcmp(get(handles.figure1,'selectiontype'),'open') && currentRoiGroup ~=0 && roinum  == size(roiData.groups(currentRoiGroup).pos,1)+1 
    
    newpos = str2num(get(handles.Position,'String'));
    roigr = roiData.groups(currentRoiGroup);
    code = roiData.groups(currentRoiGroup).code;
    if ~isempty(newpos) && ~strcmp(roigr.type,'grid')
        
        newroi = makeBinData(cat(1,roigr.pos,newpos),'label',label,'type',type,...
                'trials',roigr.activeTrials,'code',code);
        roiData.groups(currentRoiGroup) = newroi.groups;    
        setappdata(parent,'CurrentRoi',size(newroi.groups.pos,1))
        setappdata(parent,'roiData',roiData)
    end
end

SetFields(hObject, eventdata, handles)
UpdateRois(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function roiList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeActive(hObject, eventdata, handles)
% Makes the current control active in the main window

parent = getappdata(handles.figure1,'parent');

setappdata(parent,'activeControl','roiManager');




% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

if ishandle(parent)
    ph = getappdata(parent,'patchHandles');
    delete(ph(ishandle(ph)))

    setappdata(parent,'activeControl','main');
end
%-------------------------------

function    SetFields(hObject, eventdata, handles, varargin)

if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
else
    parent = varargin{1};
end

roiData = getappdata(parent,'roiData');

currentRoiGroup = getappdata(parent,'CurrentRoiGroup');
currentRoi = getappdata(parent,'CurrentRoi');

set(handles.roiNumber,'string',num2str(currentRoi));
set(handles.roiList,'value',currentRoi+1);



if ishandle(handles.figure1)
    if isempty(currentRoi) || currentRoi == 0 
        set(handles.Position,'String','[ ]')
    else
        set(handles.Position,'String',['[ ',num2str(roiData.groups(currentRoiGroup).pos(currentRoi,:)),' ]'])
    end    

%     bgstr  = {roiData.groups.label};
%     set(handles.roiGroupList,'Value',currentRoiGroup+1);
%     set(handles.roiGroupList,'String',cat(2,{'New Roi Group'},bgstr));
end

if currentRoiGroup ~=0 && ishandle(handles.figure1)
    
    %update the roi manager windo if it exists
    
	label = roiData.groups(currentRoiGroup).label;
    
    roiType = roiData.groups(currentRoiGroup).type;
    roigr = roiData.groups(currentRoiGroup);
    
        set(handles.label,'String',label)

    if currentRoi > size(roigr.pos,1)
        currentRoi  = 0;
        setappdata(parent,'CurrentRoi',currentRoi)
    end

    stringpos = mat2cell(num2str(roigr.pos,2),ones(size(roigr.pos,1),1));

    if ~isempty(stringpos) && ~strcmp(roiData.groups(currentRoiGroup).type,'poly')
        roinums = strcat(mat2cell(num2str(roigr.binnums(:)),ones(size(roigr.pos,1),1)),{'.  '});
    else
        roinums = {};
    end
    
    stringpos = strcat(roinums,stringpos);

    switch roiType

        case 'grid'

            colLabel = '#    xlim1    xlim2    ylim1    ylim2';
            set(handles.type,'Value',4)

            set(handles.label,'String',roigr.label)
        case 'rect'
            colLabel = '#    xlim1    xlim2    ylim1    ylim2';
            set(handles.type,'Value',1)  

            set(handles.label,'String',roigr.label)
        case 'ellipse'
            colLabel = '#    X  Y  r1  r2  theta';
            set(handles.type,'Value',2)
            set(handles.label,'String',roigr.label)
        case 'poly'
            set(handles.type,'Value',3)
            colLabel = 'Vertex     X  	Y ';
        case 'undefined'
            colLabel = '  ';

%             return
        otherwise
            error('Unrecognized roi type: %s',roiType)
    end
    if all(cellfun(@isempty,stringpos))
        stringpos={};
    end
    if strcmp(roiData.groups(currentRoiGroup).type,'poly')
        set(handles.roiList,'String',cat(1,colLabel,stringpos,{'New Vertex'}));
    else
        set(handles.roiList,'String',cat(1,colLabel,stringpos,{'New Roi'}));
    end
elseif ishandle(handles.figure1)
    set(handles.label,'string','none selected')
    set(handles.Position,'string','[ ]')
%     set(handles.Position,'string','[ ]')
end



function Position_Callback(hObject, eventdata, handles)
% hObject    handle to Position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Position as text
%        str2double(get(hObject,'String')) returns contents of Position as a double

parent = getappdata(handles.figure1,'parent');
currRoi = getappdata(parent,'CurrentRoi');
currentRoiGroup  = getappdata(parent,'CurrentRoiGroup');

if currRoi == 0|| currentRoiGroup ==0
    return
end

roiData = getappdata(parent,'roiData');

newpos = str2num(get(handles.Position,'String'));

% if strcmp(roiData.groups(currenTroiGroup),'grid')
    

if size(newpos,2) == size(roiData.groups(currentRoiGroup).pos,2) ||size(roiData.groups(currentRoiGroup).pos,2) == 0
    roiData.groups(currentRoiGroup).pos(currRoi,:)= newpos;
    setappdata(parent,'roiData',roiData)
else    
    set(handles.Position,'String', num2str(   roiData.groups(currentRoiGroup).pos(currRoi,:) ));    
end

if currRoi <= size(roiData.groups(currentRoiGroup).pos,1)
     UpdateRois(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function Position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roiNumber_Callback(hObject, eventdata, handles)
% hObject    handle to roiNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roiNumber as text
%        str2double(get(hObject,'String')) returns contents of roiNumber as a double


parent = getappdata(handles.figure1,'parent');
currentRoiGroup = getappdata(parent,'CurrentRoiGroup');
roiData = getappdata(parent,'roiData');

if currentRoiGroup == 0
    return
end
roinum = str2double(get(handles.roiNumber,'String'));
if ~isempty(roinum)  && roinum <= size(roiData.groups(currentRoiGroup).pos,1)

    setappdata(parent,'CurrentRoi',roinum)
else
    currroi = getappdata(parent,'CurrentRoi');
    set(handles.roiNumber,'String',num2str(currroi)); 
        
end

UpdateRois(hObject, eventdata, handles)


%-----------------

function UpdateRois(hObject, eventdata, handles,varargin)

if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
else
    parent = varargin{1};
end

currentRoiGroup = getappdata(parent,'CurrentRoiGroup');
roiData = getappdata(parent,'roiData');
phandles = guidata(parent);
cla(phandles.axes2)

while currentRoiGroup > length(roiData.groups)
    currentRoiGroup  = currentRoiGroup - 1;
    setappdata(parent,'CurrentRoiGroup',currentRoiGroup )    
end

if currentRoiGroup ~=0
    code = roiData.groups(currentRoiGroup).code;
    if strcmp(roiData.groups(currentRoiGroup).type,'grid')
        rect = roiData.groups(currentRoiGroup).pos;
        
        newroigr = makeBinData(rect2grid(rect),'type','grid','code',code,...
            'label',roiData.groups(currentRoiGroup).label,'trials',roiData.groups(currentRoiGroup).activeTrials);
%         newroigr.groups.type = 'grid';
    else
        newroigr = makeBinData(roiData.groups(currentRoiGroup).pos,'type',roiData.groups(currentRoiGroup).type,'code',code,...
            'label',roiData.groups(currentRoiGroup).label,'trials',roiData.groups(currentRoiGroup).activeTrials);
    end    
    roiData.groups(currentRoiGroup) = newroigr.groups;
    setappdata(parent,'roiData',roiData);
end

ClearPatch(hObject, eventdata, handles,varargin{:});
SetFields(hObject, eventdata, handles,varargin{:});
DrawPatch(hObject, eventdata, handles,varargin{:});
if hObject == handles.figure1
    figure(handles.figure1)
end

%Now Update trials now

trfun = getappdata(parent,'trialManagerFunctions');
if ~isempty(trfun)    
    trfun.updateFixations(parent);
end
    
% --- Executes during object creation, after setting all properties.
function roiNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function ph = DrawPatch(hObject, eventdata, handles,varargin)
% hObject    handle to roiNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

defaultAlpha = .5;

if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
    currentRoiGroups = get(handles.roiGroupList,'value')-1;
else
    parent = varargin{1};
    if length(varargin) > 1
        currentRoiGroups = varargin{2};
    else
        currentRoiGroups = 0;
    end
end

if nargin > 3 && ishandle(varargin{1})
    varargin(1) = [];
end

% currentRoiGroup = getappdata(parent,'CurrentRoiGroup');
% currentRoi = getappdata(parent,'CurrentRoi');
if currentRoiGroups(1) == 0
    currentRoiGroups = 0;
    getRois = 0;
else
    if ishandle(handles.figure1)
        getRois = get(handles.roiList,'value')-1;
    else
        getRois = 0;
    end
    if ~isequal(getRois,0)
        getRois( getRois == 0 | getRois >= length(get(handles.roiList,'string'))-1 ) = [];    
    end
end

roiData = getappdata(parent,'roiData');

if (length(roiData.groups) == 1 && isempty(roiData.groups(1).pos) )||...
        ( ~isequal(currentRoiGroups,0) && all(cellfun(@isempty, {roiData.groups(currentRoiGroups).pos})))
    return
end
% screenData = getappdata(parent,'screenData');
% screenres = screenData.res;

if isequal(currentRoiGroups,0)
    currentRoiGroups = 1:length(roiData.groups);
% else
%     currentRoiGroups = currentRoiGroup;
end
phandles = guidata(parent);
hold(phandles.axes2,'on')
for currentRoiGroup = currentRoiGroups 
    roigr = roiData.groups(currentRoiGroup);

%     type = roigr.type;
    pos = roigr.pos;
    if isempty(pos)
        continue
    end
%     if currentRoi == 0
%         getRoi = 1:size(pos,1);
%     else
%         getRoi = currentRoi;
%     end
    if isequal(getRois , 0 )
        getRoi = 1:size(pos,1);
%     else
%         getRoi = currentRoi;
    else
        getRoi = getRois;
    end

    if isempty(varargin) || ~isnumeric(varargin{1})
        C = (1:length(getRoi))./length(getRoi);
    end
    phold = getappdata(parent,'patchHandles');
    
%     axes(phandles.axes2)
    ph = cat(2,phold,roigr.patch(roigr,getRoi,C,'FaceAlpha',defaultAlpha,'parent',phandles.axes2)); 
%     ph = cat(2,phold,patch(vx*screenres(1),vy*screenres(2),C,'FaceAlpha',defaultAlpha,'parent',phandles.axes2));
%     axis(phandles.axes2,[0 screenres(1) 0 screenres(2)])
    setappdata(parent,'patchHandles',ph);
end

hold(phandles.axes2,'off')

%----------------

function ph = ClearPatch(hObject, eventdata, handles,varargin)

if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
else
    parent = varargin{1};
end

ph = getappdata(parent,'patchHandles');

delete(ph(ishandle(ph)))

setappdata(parent,'patchHandles',[]);
       
        


% --------------------------------------------------------------------
function RoiGrListContext_Callback(hObject, eventdata, handles)
% hObject    handle to RoiGrListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RemoveRoiGr_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveRoiGr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
% currentRoiGroup = getappdata(parent,'CurrentRoiGroup');
roiData = getappdata(parent,'roiData');

selected = get(handles.roiGroupList,'value')-1;
selected(selected == 0) = [];

% if currentRoiGroup == 0 
%     return
% end


roiData.groups(selected) = [];

setappdata(parent,'roiData',roiData);
setappdata(parent,'CurrentRoi',0);
currentBinGroup = selected(1)-1;
setappdata(parent,'CurrentRoiGroup',currentBinGroup );

bgstr  = {roiData.groups.label};
set(handles.roiGroupList,'Value',currentBinGroup+1);
set(handles.roiGroupList,'String',cat(2,{'New ROI Group'},bgstr));

UpdateRois(hObject, eventdata, handles)


% --------------------------------------------------------------------
function RemoveRoi_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
currentRoiGroup = getappdata(parent,'CurrentRoiGroup');
% currentRoi = getappdata(parent,'CurrentRoi');
roiData = getappdata(parent,'roiData');

pos = roiData.groups(currentRoiGroup).pos;
type = roiData.groups(currentRoiGroup).type;
label = roiData.groups(currentRoiGroup).label;
trials = roiData.groups(currentRoiGroup).activeTrials;

if currentRoiGroup == 0 
    return
end

selected = get(handles.roiList,'value');
selected(selected == 0) = [];
pos(selected-1,:) = [];

if strcmp(type,'grid')
    pos = rect2grid(pos);
end

code = roiData.groups(currentRoiGroup).code;
newroidata = makeBinData(pos,'type',type,'label',label,'trials',trials,'code',code);

roiData.groups(currentRoiGroup) = newroidata.groups;

setappdata(parent,'roiData',roiData);
setappdata(parent,'CurrentRoi',selected(1)-2);

UpdateRois(hObject, eventdata, handles)


% --------------------------------------------------------------------
function RoiListContext_Callback(hObject, eventdata, handles)
% hObject    handle to RoiListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function importSelectedBinMenu_Callback(hObject, eventdata, handles)
% hObject    handle to importSelectedBinMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

binData = getappdata(parent,'binData');
roiData = getappdata(parent,'roiData');

bmh = getappdata(parent,'BinManager');
if ~ishandle(bmh)
    return
end

bmhandles = guidata(bmh);
selected = get(bmhandles.binGroupList,'value')-1;
selected(selected==0) = [];

binData.groups = binData.groups(selected);
roiData = makeBinData(roiData,binData);

setappdata(parent,'roiData',roiData);

UpdateRois(hObject, eventdata, handles)

% --------------------------------------------------------------------
function importMenu_Callback(hObject, eventdata, handles)
% hObject    handle to importMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RemoveBinGr_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveRoiGr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RemoveBin_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function BinGrListContext_Callback(hObject, eventdata, handles)
% hObject    handle to RoiGrListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function BinListContext_Callback(hObject, eventdata, handles)
% hObject    handle to RoiListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


MakeActive(hObject, eventdata, handles)
UpdateRois(hObject, eventdata, handles)
