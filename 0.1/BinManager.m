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

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% Edit the above text to modify the response to help BinManager

% Last Modified by GUIDE v2.5 08-Apr-2011 17:02:28

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

children = getappdata(parent,'children');
children(end+1) = handles.figure1;
setappdata(parent,'children',children)

setappdata(parent,'CurrentBin',0);
% currentBinGroup = getappdata(parent,'CurrentBinGroup');

% binManagerFunctions.update = @SetFields;


currentBinGroup = getappdata(parent,'CurrentBinGroup');
binData = getappdata(parent,'binData');

if isempty(currentBinGroup) 
    setappdata(parent,'CurrentBinGroup',0) ;
end

types = {'rect','ellipse','poly','grid','simplex'};
setappdata(handles.figure1,'types',types);

if isempty(binData)
    
    [binData,binfunctions] = makeBinData([],'label','Bin Group 1','type',types{ get( handles.type,'value' ) },'code',1 );
    setappdata(parent,'binData',binData)
    currentBinGroup = 1;
%     binData.codeincr = 0;
    setappdata(parent,'CurrentBinGroup',currentBinGroup)
    
end
[qq,binfunctions] = makeBinData([]);

binManagerFunctions = binfunctions;
binManagerFunctions.draw = @(varargin)DrawPatch([],[],handles,varargin{:});
binManagerFunctions.clear = @(varargin)ClearPatch([],[],handles,varargin{:});
binManagerFunctions.update = @(varargin)UpdateBins([],[],handles,varargin{:});

setappdata(parent,'binManagerFunctions',binManagerFunctions);

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

if length(binval) > 1 %&&  get(handles.displayCheck,'value')

%     ClearPatch(hObject, eventdata, handles);
%     DrawPatch(hObject, eventdata, handles);
    UpdateBins(hObject, eventdata, handles)

    return
end

setappdata(parent,'CurrentBin',0);
set(handles.binList,'value',1);
currentBinGroup = binval-1;
setappdata(parent,'CurrentBinGroup',currentBinGroup)    
if strcmp(get(handles.figure1,'selectiontype'),'open') && binval == 1  
    
    currBinGr = length(binData.groups)+1;
    
    setappdata(parent,'CurrentBinGroup',currBinGr);
%     setappdata(parent,'CurrentBin',0);
    
    types = getappdata(handles.figure1,'types');
    type = types{get(handles.type,'Value')};
    binData = makeBinData(binData,[],'type',type,'label',sprintf('Bin Group %i',currBinGr));
    
%     binData.groups(currBinGr) = newbin.groups;
    setappdata(parent,'binData',binData)
    
end

bgstr  = {binData.groups.label};
set(handles.binGroupList,'Value',currentBinGroup+1);
set(handles.binGroupList,'String',cat(2,{'New Bin Group'},bgstr));

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

% if length(binnum) > 1 && get(handles.displayCheck,'value')
%     ClearPatch(hObject, eventdata, handles);
%     DrawPatch(hObject, eventdata, handles);
%     return
% end

if currentBinGroup ~=0 && ~isempty(binnum(1)) && binnum(1)  <= binData.groups(currentBinGroup).nbin
    setappdata(parent,'CurrentBin',binnum(1))
    UpdateBins(hObject, eventdata, handles,[],false)
elseif strcmp(get(handles.figure1,'selectiontype'),'open') && currentBinGroup ~=0 && binnum(1)  ==  binData.groups(currentBinGroup).nbin+1 
    
    newpos = str2num(get(handles.Position,'String'));
    bingr = binData.groups(currentBinGroup);
    code = binData.groups(currentBinGroup).code;
    if ~isempty(newpos) && ~strcmp(bingr.type,'grid')
        
        newbin = makeBinData(cat(1,bingr.pos,newpos),'label',label,'type',type,...
                'trials',bingr.activeTrials,'code',code);
        binData.groups(currentBinGroup) = newbin.groups;    
        setappdata(parent,'CurrentBin',size(newbin.groups.pos,1))
        setappdata(parent,'binData',binData)
    end
    UpdateBins(hObject, eventdata, handles)
end


SetFields(hObject, eventdata, handles)
% UpdateBins(hObject, eventdata, handles)

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

%%% Sets the listbox and text window displays for the current bin group 

if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
else
    parent = varargin{1};
end

binData = getappdata(parent,'binData');

currentBinGroup = getappdata(parent,'CurrentBinGroup');
currentBin = getappdata(parent,'CurrentBin');

set(handles.binNumber,'string',num2str(currentBin));
% set(handles.binList,'value',currentBin+1);

% bgstr  = {binData.groups.label};



if ishandle(handles.figure1)
    if isempty(currentBin) || currentBin == 0 || isempty(binData.groups(currentBinGroup).pos)  
        set(handles.Position,'String','[ ]')
    else
        switch binData.groups(currentBinGroup).type
            case {'tri','simplex'}
            set(handles.Position,'String',['[ ',num2str(binData.groups(currentBinGroup).pos.tri(currentBin,:)),' ]'])
                
            otherwise
            set(handles.Position,'String',['[ ',num2str(binData.groups(currentBinGroup).pos(currentBin,:)),' ]'])
        end
    end    

% bgstr  = {binData.groups.label};
%     set(handles.binGroupList,'Value',currentBinGroup+1);
%     set(handles.binGroupList,'String',cat(2,{'New Bin Group'},bgstr));
end

if currentBinGroup ~=0 && ishandle(handles.figure1)
    
    %update the bin manager window if it exists
    
	label = binData.groups(currentBinGroup).label;
    
    binType = binData.groups(currentBinGroup).type;
    bingr = binData.groups(currentBinGroup);
    
    if ismember(bingr.type,{'tri','simplex'}) && ~isempty(bingr.pos)
        bingr.pos = bingr.pos.tri;              %For now only vertex numbers will be displayed.
    end
    
    
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
        case {'tri','simplex'}
            set(handles.type,'Value',5)
            colLabel = 'Vertices  ';
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
    
    val = get(handles.binList,'value');
    if val > length(cat(1,colLabel,stringpos,{'New Bin'}))
        set(handles.binList,'value',length(cat(1,colLabel,stringpos,{'New Bin'})));
    end
    
elseif ishandle(handles.figure1)
    set(handles.label,'string','none selected')
    set(handles.Position,'string','[ ]')
    set(handles.binList,'string','')
    set(handles.binList,'value',0)
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
if ismember(binData.groups(currentBinGroup),{'tri','simplex'})
    return
end


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
if ~isempty(binnum)  && binnum <=  binData.groups(currentBinGroup).nbin

    setappdata(parent,'CurrentBin',binnum)
else
    currbin = getappdata(parent,'CurrentBin');
    set(handles.binNumber,'String',num2str(currbin)); 
        
end

set(handles.binGroupList,'Value',currentBinGroup+1);

UpdateBins(hObject, eventdata, handles)


%-----------------

function UpdateBins(hObject, eventdata, handles,varargin)



% full_update = false;

if ishandle(handles.figure1)  %Get parent figure handles
    parent = getappdata(handles.figure1,'parent');
else
    parent = varargin{1};
end

if nargin < 5
    full_update = true;
else
    full_update = false;
    varargin(1:2) = [];
end

currentBinGroup = getappdata(parent,'CurrentBinGroup');

currentBinGroup = currentBinGroup(1); % Handle case when multiple groups are active 

binData = getappdata(parent,'binData');
phandles = guidata(parent);
% cla(phandles.axes2)

while currentBinGroup > length(binData.groups)
    currentBinGroup  = currentBinGroup - 1;
    setappdata(parent,'CurrentBinGroup',currentBinGroup )    
end

if full_update
   

    if currentBinGroup ~=0
        code = binData.groups(currentBinGroup).code;
        if strcmp(binData.groups(currentBinGroup).type,'grid')
            rect = binData.groups(currentBinGroup).pos;

            newbingr = makeBinData(rect2grid(rect),'type','grid','code',code,...
                'label',binData.groups(currentBinGroup).label,'trials',binData.groups(currentBinGroup).activeTrials);
    %         newbingr.groups.type = 'grid';
        else
            newbingr = makeBinData(binData.groups(currentBinGroup).pos,'type',binData.groups(currentBinGroup).type,'code',code,...
                'label',binData.groups(currentBinGroup).label,'trials',binData.groups(currentBinGroup).activeTrials);
        end    
        binData.groups(currentBinGroup) = newbingr.groups;
        setappdata(parent,'binData',binData);
    end

    %Now Update trials 

    trfun = getappdata(parent,'trialManagerFunctions');
    if ~isempty(trfun)    
        trfun.updateFixations(parent);
    end
end
ClearPatch(hObject, eventdata, handles,varargin{:});

if get(handles.displayCheck,'value')
    SetFields(hObject, eventdata, handles,varargin{:});
    DrawPatch(hObject, eventdata, handles,varargin{:});
    if hObject == handles.figure1
        activefigure(handles.figure1)
    end
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
defaultShading = 'flat';


if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
    currentBinGroups = get(handles.binGroupList,'value')-1;
elseif ishandle(varargin{1})
    parent = varargin{1};
    if length(varargin) < 2
        currentBinGroups = 0;
    end
else 
    parent = gcf;
end
 
 
if nargin > 3 && ~isempty(varargin{1}) && ishandle(varargin{1})
    varargin(1) = [];
end

 if length(varargin) > 1
        currentBinGroups = varargin{2};
 end

if length(varargin)> 1 && ~iscell(varargin{1}) && all(~ishandle(varargin{1}))
    varargin{1} = {varargin{1}};
end

% currentBinGroup = getappdata(parent,'CurrentBinGroup');
% currentBin = getappdata(parent,'CurrentBin');
if currentBinGroups(1) == 0
    currentBinGroups = 0;
    getBins = 0;
else
    if length(varargin) > 2
        getBins = varargin{3};
    elseif ishandle(handles.figure1)       
        getBins = get(handles.binList,'value')-1;
    else
        getBins = 0;
    end
    if ~isequal(getBins,0)
        getBins( getBins == 0 | getBins >= length(get(handles.binList,'string'))-1 ) = [];    
    end
end

binData = getappdata(parent,'binData');

if (length(binData.groups) == 1 &&   ~isequal(currentBinGroups,0) && binData.groups(currentBinGroups).nbin==0 )||...
        ( ~isequal(currentBinGroups,0) && all(cellfun(@isempty, {binData.groups.pos})))
    fprintf('\nSelected bin groups are empty...')
    return
end
% screenData = getappdata(parent,'screenData');
% screenres = screenData.res;

if isequal(currentBinGroups,0)
    currentBinGroups = 1:length(binData.groups);
end


if isempty(currentBinGroups)
    fprintf('\nNo active bin groups to plot...')
    return
end

phandles = guidata(parent);
hold(phandles.axes2,'on')

plotincr = 0;   
for i = 1:length(currentBinGroups) 
    bingr = binData.groups(currentBinGroups(i));

%     type = bingr.type;
    pos = bingr.pos;
    if isempty(pos)
        continue
    end
%     if currentBin == 0
%         getBin = 1:size(pos,1);
%     else
%         getBin = currentBin;
%     end
    if isequal(getBins , 0 )
        if ismember(bingr.type,{'tri','simplex'})
            getBin = 1:size(pos.tri,1);
        elseif ismember(bingr.type,{'poly'})
            getBin = 1;
        else
            getBin = 1:size(pos,1);
        end
%             else
%         getBin = currentBin;
    else
        getBin = getBins;
    end
    if isempty(getBin)
        continue
    end
    if isempty(varargin) 
        C = ((1:length(getBin)) + plotincr)./length(getBin);
        crange(i,:) = [min(C), max(C)];
    else
        C = varargin{1}{i}(:)';
        ccat = cat(1,varargin{1}{:});
        crange = [min(ccat), max(ccat)];
    end
    phold = getappdata(parent,'patchHandles');
    
%     axes(phandles.axes2)
    ph = cat(2,phold,bingr.patch(bingr,getBin,C,'FaceAlpha',defaultAlpha,'parent',phandles.axes2));
%     crange = [min(C) max(C)];
    if diff(crange) < 1e-10
        crange = [0 1];
    end
%     ph = cat(2,phold,patch(vx*screenres(1),vy*screenres(2),C,'FaceAlpha',defaultAlpha,'parent',phandles.axes2));
%     axis(phandles.axes2,[0 screenres(1) 0 screenres(2)])
    setappdata(parent,'patchHandles',ph);
    
    plotincr = plotincr+length(getBin);
end

caxis(phandles.axes2, [min(crange(:,1)) max(crange(:,2))])

shading(phandles.axes2,defaultShading)

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
% currentBinGroup = getappdata(parent,'CurrentBinGroup');
binData = getappdata(parent,'binData');

selected = get(handles.binGroupList,'value')-1;
selected(selected == 0) = [];

% if currentBinGroup == 0 
%     return
% end


binData.groups(selected) = [];

setappdata(parent,'binData',binData);
setappdata(parent,'CurrentBin',0);
currentBinGroup = selected(1)-1;
setappdata(parent,'CurrentBinGroup',currentBinGroup );

bgstr  = {binData.groups.label};
set(handles.binGroupList,'Value',currentBinGroup+1);
set(handles.binGroupList,'String',cat(2,{'New Bin Group'},bgstr));

UpdateBins(hObject, eventdata, handles)


% --------------------------------------------------------------------
function RemoveBin_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
currentBinGroup = getappdata(parent,'CurrentBinGroup');
% currentBin = getappdata(parent,'CurrentBin');
binData = getappdata(parent,'binData');

pos = binData.groups(currentBinGroup).pos;
type = binData.groups(currentBinGroup).type;
label = binData.groups(currentBinGroup).label;
trials = binData.groups(currentBinGroup).activeTrials;

if currentBinGroup == 0 
    return
end

selected = get(handles.binList,'value');
selected(selected == 0) = [];

if ismember(type,{'tri','simplex'})
    pos.tri(selected-1,:) = [];
elseif strcmp(type,'grid')    
    pos = rect2grid(pos);
else
    pos(selected-1,:) = [];

end

code = binData.groups(currentBinGroup).code;
newbindata = makeBinData(pos,'type',type,'label',label,'trials',trials,'code',code);

binData.groups(currentBinGroup) = newbindata.groups;

setappdata(parent,'binData',binData);
setappdata(parent,'CurrentBin',selected(1)-2);

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


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
UpdateBins(hObject, eventdata, handles,[],false)


% --- Executes on button press in displayCheck.
function displayCheck_Callback(hObject, eventdata, handles)
% hObject    handle to displayCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayCheck

parent = getappdata(handles.figure1,'parent');
rmh = getappdata(parent,'RegManager');
if ~isempty(rmh) && ishandle(rmh) && get(handles.displayCheck,'value')
    rmhandles = guidata(rmh);
    set(rmhandles.displayCheck,'value',0)
end

UpdateBins(hObject, eventdata, handles)


% --- Executes on button press in edit_bins_check.
function edit_bins_check_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bins_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_bins_check


% --------------------------------------------------------------------
function copybingr_Callback(hObject, eventdata, handles)
% hObject    handle to copybingr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



parent = getappdata(handles.figure1,'parent');
currentBinGroup = getappdata(parent,'CurrentBinGroup');
% currentBin = getappdata(parent,'CurrentBin');
binData = getappdata(parent,'binData');
vals = get(handles.binGroupList,'value');
vals(vals==1) = [];
if ~isempty(vals)
    for i = 1:length(vals)
        binData.groups(end+1) = binData.groups(vals(i)-1);
        binData.groups(end).code = binData.codeincr+1;
        binData.codeincr = binData.codeincr+1;
    end
    setappdata(parent,'binData',binData);
end
