function varargout = BinManager(varargin)
% BINMANAGER M-file for BinManager.fig
%      BINMANAGER, by itself, creates a new BINMANAGER or raises the existing
%      singleton*.
%
%      H = BINMANAGER returns the handle to a new BINMANAGER or the handle to
%      the existing singleton*.
%
%      BINMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BINMANAGER.M with the given input
%      arguments.
%
%      BINMANAGER('Property','Value',...) creates a new BINMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BinManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BinManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BinManager

% Last Modified by GUIDE v2.5 19-Dec-2007 14:42:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BinManager_OpeningFcn, ...
                   'gui_OutputFcn',  @BinManager_OutputFcn, ...
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


% --- Executes just before BinManager is made visible.
function BinManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BinManager (see VARARGIN)

% Choose default command line output for BinManager
handles.output = hObject;

parent =varargin{1};
setappdata(handles.figure1,'parent',parent );
MakeActive(hObject, eventdata, handles)

setappdata(parent,'CurrentBin',0);
% currentBinGroup = getappdata(parent,'CurrentBinGroup');

binManagerFunctions.draw = @DrawPatch;
binManagerFunctions.clear = @ClearPatch;
binManagerFunctions.update = @UpdateBins;
% binManagerFunctions.update = @SetFields;

setappdata(parent,'binManagerFunctions',binManagerFunctions);

currentBinGroup = getappdata(parent,'CurrentBinGroup');
binData = getappdata(parent,'binData');

if isempty(currentBinGroup)
    setappdata(parent,'CurrentBinGroup',0) ;
end

types = {'rect','ellipse','poly','grid'};
setappdata(handles.figure1,'types',types);

if isempty(binData)
    
    binData = makeBinData([],'label','Bin Group 1','type',types{ get( handles.type,'value' ) } );
    setappdata(parent,'binData',binData)
    currentBinGroup = 1;
    setappdata(parent,'CurrentBinGroup',currentBinGroup)
    
end

UpdateBins(hObject, eventdata, handles)
% if currentBinGroup ~= 0
%     SetFields(hObject, eventdata, handles);
% end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BinManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BinManager_OutputFcn(hObject, eventdata, handles) 
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
currBinGr = getappdata(parent,'CurrentBinGroup');

if currBinGr ~= 0 
    binData = getappdata(parent,'binData');
    binData.groups(currBinGr).label = get(handles.label,'String');
    setappdata(parent,'binData',binData)
end


SetFields(hObject, eventdata, handles)
UpdateBins(hObject, eventdata, handles)


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
currBinGr = getappdata(parent,'CurrentBinGroup');

binData = getappdata(parent,'binData');
val = get(handles.type,'value');

if currBinGr ~= 0 
%      strs = get(handles.type,'string');

    if ~isempty(binData.groups(currBinGr).pos)
        set(handles.type,'value',find(ismember(typeLabels,binData.groups(currBinGr).type)));
    else
        binData.groups(currBinGr).type = typeLabels{val};
        setappdata(parent,'binData',binData)
    end        
end

SetFields(hObject, eventdata, handles)
UpdateBins(hObject, eventdata, handles)


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


% --- Executes on selection change in binGroupList.
function binGroupList_Callback(hObject, eventdata, handles)
% hObject    handle to binGroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns binGroupList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from binGroupList
MakeActive(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');

binData = getappdata(parent,'binData');

binval = get(handles.binGroupList,'Value');


setappdata(parent,'CurrentBin',0);
setappdata(parent,'CurrentBinGroup',binval-1)    
if binval == 1 && strcmp(get(handles.figure1,'selectiontype'),'open')
    
    currBinGr = length(binData.groups)+1;
    
    setappdata(parent,'CurrentBinGroup',currBinGr);
%     setappdata(parent,'CurrentBin',0);
    
    types = getappdata(handles.figure1,'types');
    type = types{get(handles.type,'Value')};
    newbin = makeBinData([],'type',type,'label',sprintf('Bin Group %i',currBinGr));
    
    binData.groups(currBinGr) = newbin.groups;
    setappdata(parent,'binData',binData)
    
end
    
SetFields(hObject, eventdata, handles)
UpdateBins(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function binGroupList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binGroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in binList.
function binList_Callback(hObject, eventdata, handles)
% hObject    handle to binList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns binList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from binList


parent = getappdata(handles.figure1,'parent');
binData = getappdata(parent,'binData');

currentBinGroup = getappdata(parent,'CurrentBinGroup');
binnum = get(handles.binList,'Value')-1;
types = getappdata(handles.figure1,'types');
type = types{get(handles.type,'value')};
label = get(handles.label,'string');

if currentBinGroup ~=0 && ~isempty(binnum) && binnum  <= size(binData.groups(currentBinGroup).pos,1)
    setappdata(parent,'CurrentBin',binnum)
elseif currentBinGroup ~=0 && binnum  == size(binData.groups(currentBinGroup).pos,1)+1
    
    newpos = str2num(get(handles.Position,'String'));
    bingr = binData.groups(currentBinGroup);
    
    if ~isempty(newpos) && ~strcmp(bingr.type,'grid')
        
        newbin = makeBinData(cat(1,bingr.pos,newpos),'label',label,'type',type,...
                'trials',bingr.activeTrials);
        binData.groups(currentBinGroup) = newbin.groups;    
        setappdata(parent,'CurrentBin',size(newbin.groups.pos,1))
        setappdata(parent,'binData',binData)
    end
end

SetFields(hObject, eventdata, handles)
UpdateBins(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function binList_CreateFcn(hObject, eventdata, handles)
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

setappdata(parent,'activeControl','binManager');


% function updatePlot(hObject,eventdata,roiInstance,pivot,btnDownPt,axh)
% function updatePlot(hObject,eventdata,handles)
%    
% plottype = 'g'; %arguments to plot for drawing bin lines
% 
% % Simply calls the plotting function handle associated with this bin group
% 
% parent = getappdata(handles.figure1,'parent');
% phandles = guidata(parent);
% axes(phandles.axes2)
% gridData = getappdata(parent,'gridData');
% % currim = getappdata(parent,'CurrentImage');
% % trialData = getappdata(parent,'trialData');
% % screenData = getappdata(parent,'screenData');
% 
% 
% currentgr = getappdata(parent,'CurrentBinGroup');
% 
% phandles = guidata(parent); %handles of the parent GUI
% 
% 
% plothandles = getappdata(handles.figure1,'plotHandles');
% 
% delete(plothandles(ishandle(plothandles(:))));
% 
% axes(phandles.axes2)
% 
% axis(phandles.axes2,axis(phandles.axes1));
% ax1 = get(phandles.axes1);
% set(phandles.axes2,'units',ax1.Units,'position',ax1.Position,'Ydir',ax1.YDir)
% 
% hold off
% plothandles = gridData.groups(currentgr).plot(gridData.groups(currentgr),plottype);
% % hold off
% axis off
% % axes(phandles.axes1)
% 
% setappdata(handles.figure1,'plotHandles',plothandles)




% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ph = getappdata(hObject,'patchHandles');
delete(ph(ishandle(ph)))

%-------------------------------

function    SetFields(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');

binData = getappdata(parent,'binData');

currentBinGroup = getappdata(parent,'CurrentBinGroup');
currentBin = getappdata(parent,'CurrentBin');

set(handles.binNumber,'string',num2str(currentBin));
set(handles.binList,'value',currentBin+1);

bgstr  = {binData.groups.label};


if isempty(currentBin) || currentBin == 0 
    set(handles.Position,'String','[ ]')
else
    set(handles.Position,'String',['[ ',num2str(binData.groups(currentBinGroup).pos(currentBin,:)),' ]'])
end    

set(handles.binGroupList,'Value',currentBinGroup+1);
set(handles.binGroupList,'String',cat(2,{'New Bin Group'},bgstr));
if currentBinGroup ~=0

	label = binData.groups(currentBinGroup).label;
    
    binType = binData.groups(currentBinGroup).type;
    bingr = binData.groups(currentBinGroup);
    
    set(handles.label,'String',label)
    

    if currentBin > size(bingr.pos,1)
        currentBin  = 0;
        setappdata(parent,'CurrentBin',currentBin)
    end

    stringpos = mat2cell(num2str(bingr.pos,2),ones(size(bingr.pos,1),1));

    if ~isempty(stringpos) && ~strcmp(binData.groups(currentBinGroup).type,'poly')
        binnums = strcat(mat2cell(num2str(bingr.binnums(:)),ones(size(bingr.pos,1),1)),{'.  '});
    else
        binnums = {};
    end
    
    stringpos = strcat(binnums,stringpos);

    switch binType

        case 'grid'

            colLabel = '#    xlim1    xlim2    ylim1    ylim2';
            set(handles.type,'Value',4)

            set(handles.label,'String',bingr.label)
        case 'rect'
            colLabel = '#    xlim1    xlim2    ylim1    ylim2';
            set(handles.type,'Value',1)  

            set(handles.label,'String',bingr.label)
        case 'ellipse'
            colLabel = '#    X  Y  r1  r2  theta';
            set(handles.type,'Value',2)
            set(handles.label,'String',bingr.label)
        case 'poly'
            set(handles.type,'Value',3)
            colLabel = 'Vertex     X  	Y ';
        case 'undefined'
            colLabel = '  ';

%             return
        otherwise
            error('Unrecognized bin type: %s',binType)
    end
    if all(cellfun(@isempty,stringpos))
        stringpos={};
    end
    if strcmp(binData.groups(currentBinGroup).type,'poly')
        set(handles.binList,'String',cat(1,colLabel,stringpos,{'New Vertex'}));
    else
        set(handles.binList,'String',cat(1,colLabel,stringpos,{'New Bin'}));
    end
else
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
currBin = getappdata(parent,'CurrentBin');
currentBinGroup  = getappdata(parent,'CurrentBinGroup');

if currBin == 0|| currentBinGroup ==0
    return
end

binData = getappdata(parent,'binData');

newpos = str2num(get(handles.Position,'String'));

% if strcmp(binData.groups(currenTbinGroup),'grid')
    

if size(newpos,2) == size(binData.groups(currentBinGroup).pos,2) ||size(binData.groups(currentBinGroup).pos,2) == 0
    binData.groups(currentBinGroup).pos(currBin,:)= newpos;
    setappdata(parent,'binData',binData)
else    
    set(handles.Position,'String', num2str(   binData.groups(currentBinGroup).pos(currBin,:) ));    
end

if currBin <= size(binData.groups(currentBinGroup).pos,1)
     UpdateBins(hObject, eventdata, handles)
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



function binNumber_Callback(hObject, eventdata, handles)
% hObject    handle to binNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binNumber as text
%        str2double(get(hObject,'String')) returns contents of binNumber as a double


parent = getappdata(handles.figure1,'parent');
currentBinGroup = getappdata(parent,'CurrentBinGroup');
binData = getappdata(parent,'binData');

if currentBinGroup == 0
    return
end
binnum = str2double(get(handles.binNumber,'String'));
if ~isempty(binnum)  && binnum <= size(binData.groups(currentBinGroup).pos,1)

    setappdata(parent,'CurrentBin',binnum)
else
    currbin = getappdata(parent,'CurrentBin');
    set(handles.binNumber,'String',num2str(currbin)); 
        
end

UpdateBins(hObject, eventdata, handles)


%-----------------

function UpdateBins(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');
currentBinGroup = getappdata(parent,'CurrentBinGroup');
binData = getappdata(parent,'binData');
phandles = guidata(parent);
cla(phandles.axes2)

while currentBinGroup > length(binData.groups)
    currentBinGroup  = currentBinGroup - 1;
    setappdata(parent,'CurrentBinGroup',currentBinGroup )    
end

if currentBinGroup ~=0
    if strcmp(binData.groups(currentBinGroup).type,'grid')
        rect =binData.groups(currentBinGroup).pos;
        
        newbingr = makeBinData(rect2grid(rect),'type','grid',...
            'label',binData.groups(currentBinGroup).label,'trials',binData.groups(currentBinGroup).activeTrials);
%         newbingr.groups.type = 'grid';
    else
        newbingr = makeBinData(binData.groups(currentBinGroup).pos,'type',binData.groups(currentBinGroup).type,...
            'label',binData.groups(currentBinGroup).label,'trials',binData.groups(currentBinGroup).activeTrials);
    end    
    binData.groups(currentBinGroup) = newbingr.groups;
    setappdata(parent,'binData',binData);
end

ClearPatch(hObject, eventdata, handles);
SetFields(hObject, eventdata, handles);
DrawPatch(hObject, eventdata, handles);
if hObject == handles.figure1
    figure(handles.figure1)
end

% --- Executes during object creation, after setting all properties.
function binNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function ph = DrawPatch(hObject, eventdata, handles,varargin)
% hObject    handle to binNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

defaultAlpha = .5;

parent = getappdata(handles.figure1,'parent');

currentBinGroup = getappdata(parent,'CurrentBinGroup');
currentBin = getappdata(parent,'CurrentBin');

binData = getappdata(parent,'binData');

if (length(binData.groups) == 1 && isempty(binData.groups(1).pos) )||...
        ( currentBinGroup ~= 0 && isempty(binData.groups(currentBinGroup).pos))
    return
end
screenData = getappdata(parent,'screenData');
screenres = screenData.res;

if currentBinGroup == 0
    currentBinGroups = 1:length(binData.groups);
else
    currentBinGroups = currentBinGroup;
end
phandles = guidata(parent);
hold(phandles.axes2,'on')
for currentBinGroup = currentBinGroups 
    bingr = binData.groups(currentBinGroup);

    type = bingr.type;
    pos = bingr.pos;
    if isempty(pos)
        continue
    end
    if currentBin == 0
        getBin = 1:size(pos,1);
    else
        getBin = currentBin;
    end
%     switch type
% 
%         case 'grid'    
% 
%             vx = cat(2,pos(getBin ,[1,2]),pos(getBin ,[2,1]))';
% 
%             vy = pos(getBin ,[3,3,4,4])';
% 
%         case 'rect'        
% 
%             vx = cat(2,pos(getBin ,[1,2]),pos(getBin ,[2,1]))';
%             
%             vy = pos(getBin ,[3,3,4,4])';
% 
%         case 'ellipse'
%             pos = pos(getBin,:);
%             if size(pos,2) == 3
%                pos = pos(:,[1 2 3 3 ]); %Done so that the function won't have to be modified for elliptical bins
%                pos(:,end+1) = 0;
%             end
%             
%             rotangle = 0;
% %             E = [cos(rotangle), sin(rotangle); -sin(rotangle), cos(rotangle)]; %Rotation matrix   
%             T = [pos(3)*cos(rotangle), pos(4)*sin(rotangle) 0; -pos(3)*sin(rotangle), pos(4)*cos(rotangle), 0; pos(1:2),1]; %Rotation matrix   
% 
%             th = 0:.01:2*pi;
% 
% %            crc =  [kron(pos(:,3),ones(length(th),1)).*repmat(cos(th'),size(pos,1),1)...
% %                     kron(pos(:,4),ones(length(th),1)).*repmat(sin(th'),size(pos,1),1)]*E...
% %                     + kron(pos(:,1:2),ones(length(th),1));
%            crc =  [kron(pos(:,3),ones(length(th),1)).*repmat(cos(th'),size(pos,1),1)...
%                     kron(pos(:,4),ones(length(th),1)).*repmat(sin(th'),size(pos,1),1)]*E...
%                     + kron(pos(:,1:2),ones(length(th),1));
%            
%            vx = reshape(crc(:,1),length(th),size(pos,1));
%            vy = reshape(crc(:,2),length(th),size(pos,1));
% 
%         case 'poly'
% 
%             vx = pos(:,1);
%             vy = pos(:,2);
%         case 'undefined'
%             return
%         otherwise
%             error('Unrecognized bin type')
%     end
%             

    if isempty(varargin) || ~isnumeric(varargin{1})
        C = (1:length(getBin))./length(getBin);
    end
    phold = getappdata(handles.figure1,'patchHandles');
    
%     axes(phandles.axes2)
    ph = bingr.patch(bingr,getBin,C,'FaceAlpha',defaultAlpha,'parent',phandles.axes2); 
%     ph = cat(2,phold,patch(vx*screenres(1),vy*screenres(2),C,'FaceAlpha',defaultAlpha,'parent',phandles.axes2));
%     axis(phandles.axes2,[0 screenres(1) 0 screenres(2)])
    setappdata(handles.figure1,'patchHandles',ph);
end
hold(phandles.axes2,'off')

%----------------

function ph = ClearPatch(hObject, eventdata, handles,varargin)
        
ph = getappdata(handles.figure1,'patchHandles');

delete(ph(ishandle(ph)))

setappdata(handles.figure1,'patchHandles',[]);
       
        


% --------------------------------------------------------------------
function BinGrListContext_Callback(hObject, eventdata, handles)
% hObject    handle to BinGrListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RemoveBinGr_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveBinGr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
currentBinGroup = getappdata(parent,'CurrentBinGroup');
binData = getappdata(parent,'binData');

if currentBinGroup == 0 
    return
end

binData.groups(currentBinGroup) = [];

setappdata(parent,'binData',binData);
setappdata(parent,'CurrentBin',0);
setappdata(parent,'CurrentBinGroup',currentBinGroup-1);

UpdateBins(hObject, eventdata, handles)


% --------------------------------------------------------------------
function RemoveBin_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
currentBinGroup = getappdata(parent,'CurrentBinGroup');
currentBin = getappdata(parent,'CurrentBin');
binData = getappdata(parent,'binData');

pos = binData.groups(currentBinGroup).pos;
type = binData.groups(currentBinGroup).type;
label = binData.groups(currentBinGroup).label;
trials = binData.groups(currentBinGroup).activeTrials;

if currentBinGroup == 0 || currentBin == 0 
    return
end

pos(currentBin,:) = [];

if strcmp(type,'grid')
    pos = rect2grid(pos);
end

newbindata = makeBinData(pos,'type',type,'label',label,'trials',trials);

binData.groups(currentBinGroup) = newbindata.groups;

setappdata(parent,'binData',binData);
setappdata(parent,'CurrentBin',currentBin-1);

UpdateBins(hObject, eventdata, handles)


% --------------------------------------------------------------------
function BinListContext_Callback(hObject, eventdata, handles)
% hObject    handle to BinListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------

% function grid = rect2grid(rect)
% 
% %Converts the input data format for rect type to that for grid
% 
% if iscell(rect) || isempty(rect)
%     grid = rect;
%     return
% end
% 
% minmax = cat(1,min(rect,[],1),max(rect,[],1));
% 
% grid{1} = minmax([1 4 5 8]);
% 
% nbins = [length(unique(rect(:,1))),length(unique(rect(:,3)))];
% grid{2} = nbins;
% 
