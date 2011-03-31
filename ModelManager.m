function varargout = ModelManager(varargin)
% MODELMANAGER M-file for ModelManager.fig
%      MODELMANAGER, by itself, creates a new MODELMANAGER or raises the existing
%      singleton*.
%
%      H = MODELMANAGER returns the handle to a new MODELMANAGER or the handle to
%      the existing singleton*.
%
%      MODELMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODELMANAGER.M with the given input
%      arguments.
%
%      MODELMANAGER('Property','Value',...) creates a new MODELMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModelManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModelManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ModelManager

% Last Modified by GUIDE v2.5 01-Feb-2008 15:58:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModelManager_OpeningFcn, ...
                   'gui_OutputFcn',  @ModelManager_OutputFcn, ...
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


% --- Executes just before ModelManager is made visible.
function ModelManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModelManager (see VARARGIN)

% Choose default command line output for ModelManager
handles.output = hObject;

parent =varargin{1};
setappdata(handles.figure1,'parent',parent );
% MakeActive(hObject, eventdata, handles)
setappdata(parent,'ModelManager',handles.figure1);

children = getappdata(parent,'children');
children(end+1) = handles.figure1;
setappdata(parent,'children',children)

setappdata(parent,'currentRegressorGroup',0);
% currentModel = getappdata(parent,'CurrentModel');


modelManagerFunctions.update = @()Update([],[],handles);

modelManagerFunctions.fit_response = @(Y,varargin) FitButton_Callback(hObject, eventdata, handles,Y,varargin{:});

setappdata(parent,'modelManagerFunctions',modelManagerFunctions);

currentModel = getappdata(parent,'CurrentModel');
modelData = getappdata(parent,'modelData');

if isempty(currentModel)
    setappdata(parent,'CurrentModel',0) ;
end


if isempty(modelData)
    
    modelData = makeModelData([],'label','Model 1' );
%     modelData.models(:) = [];
    setappdata(parent,'modelData',modelData)
    currentModel = 1;
%     modelData(currentDataSet).codeincr = 0;
    setappdata(parent,'CurrentModel',currentModel)
    
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ModelManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


Update(hObject, eventdata, handles)
% --- Outputs from this function are returned to the command line.
function varargout = ModelManager_OutputFcn(hObject, eventdata, handles) 
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
currModel = getappdata(parent,'CurrentModel');
currentDataSet = getappdata(parent,'CurrentDataSet');
if currModel ~= 0 
    modelData = getappdata(parent,'modelData');
    modelData(currentDataSet).models(currModel).label = get(handles.label,'String');
    setappdata(parent,'modelData',modelData)
end


Update(hObject, eventdata, handles)

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



% --- Executes on selection change in modelList.
function modelList_Callback(hObject, eventdata, handles)
% hObject    handle to modelList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns modelList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modelList


parent = getappdata(handles.figure1,'parent');
modelData = getappdata(parent,'modelData');

currentDataSet = getappdata(parent,'CurrentDataSet');

selected = get(handles.modelList,'value')-1;

selected(selected ==0 ) = [];

if ~isempty(selected)   
    setappdata(parent,'CurrentModel',selected(1));
elseif isequal(get(handles.figure1,'selectiontype'),'open')
    if length(modelData) < currentDataSet
        modelData(currentDataSet) = makeModelData([],'Label',sprintf('model %i',1));
    else        
        modelData(currentDataSet) = makeModelData(modelData(currentDataSet),'Label',sprintf('model %i',length(modelData(currentDataSet).models)+1));
    end
    setappdata(parent,'modelData',modelData);
    setappdata(parent,'CurrentModel',length(modelData(currentDataSet).models));
end

Update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function modelList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in regressorGroupList.
function regressorGroupList_Callback(hObject, eventdata, handles)
% hObject    handle to regressorGroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns regressorGroupList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from regressorGroupList


parent = getappdata(handles.figure1,'parent');
modelData = getappdata(parent,'modelData');
regData = getappdata(parent,'regData');
currentModel = getappdata(parent,'CurrentModel');
currentDataSet = getappdata(parent,'CurrentDataSet'); 
selected = get(handles.regressorGroupList,'value')-1;

selected(selected ==0 ) = [];

if ~isempty(selected) && ~(isempty(currentModel) || currentModel == 0)
    modelreg = find(ismember([regData(currentDataSet).regressors.code],modelData(currentDataSet).models(currentModel).regressors));
    setappdata(parent,'CurrentRegressorGroup',modelreg(selected(1)));
    setappdata(parent,'SelectedRegressorGroups',modelreg(selected));
%     setappdata(parent,'CurrentRegressorGroup',find(ismember([regData(currentDataSet).regressors.code],modelData(currentDataSet).models(currentModel).regressors(selected))));
elseif isequal(get(handles.figure1,'selectiontype'),'open') && ~(isempty(currentModel) || currentModel == 0)
     
    
%     rmh = getappdata(parent,'RegManager');
%     rmhandles = guidata(rmh);
% 
%     selreg = get(rmhandles.regressorGroupList,'value')-1;
%     selreg(selreg == 0 ) = [];
% 
%     modelData(currentDataSet).models(currentModel).regressors = sort(cat(2,modelData(currentDataSet).models(currentModel).regressors,[regData(currentDataSet).regressors(selreg).code]));
% 
%     setappdata(parent,'modelData',modelData);
    addRegressorsMenu_Callback(hObject, eventdata, handles)
    return
end
  


Update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function regressorGroupList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function RemoveModel_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
% currentModel = getappdata(parent,'CurrentModel');
modelData = getappdata(parent,'modelData');

selected = get(handles.modelList,'value')-1;
selected(selected == 0) = [];
currentDataSet = getappdata(parent,'CurrentDataSet');

% if currentModel == 0 
%     return
% end


modelData(currentDataSet).models(selected) = [];

setappdata(parent,'modelData',modelData);
setappdata(parent,'CurrentRegressorGroup',0);
currentModel = selected(1)-1;
setappdata(parent,'CurrentModel',currentModel );
% 
% bgstr  = {modelData(currentDataSet).models.label};
 set(handles.modelList,'Value',currentModel+1);
 set(handles.regressorGroupList,'Value',1);
% set(handles.modelList,'String',cat(2,{'New Model'},bgstr));

Update(hObject, eventdata, handles)


% --------------------------------------------------------------------
function RemoveRegressorGroup_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveRegressorGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

currentModel = getappdata(parent,'CurrentModel');
modelData = getappdata(parent,'modelData');
currentDataSet = getappdata(parent,'CurrentDataSet');

if currentModel == 0 
    return
end

selected = get(handles.regressorGroupList,'value')-1;
selected(selected == 0) = [];

modelData(currentDataSet).models(currentModel).regressors(selected) = [];

setappdata(parent,'modelData',modelData);
set(handles.regressorGroupList,'value',selected(1));

Update(hObject, eventdata, handles)

% --------------------------------------------------------------------
function RegressorListContext_Callback(hObject, eventdata, handles)
% hObject    handle to RegressorGroupListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% MakeActive(hObject, eventdata, handles)
% UpdateBins(hObject, eventdata, handles)


% --- Executes on button press in displayDesign.
function displayDesign_Callback(hObject, eventdata, handles)
% hObject    handle to displayDesignmodelData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayDesign




% --- Executes on button press in FitButton.
function FitButton_Callback(hObject, eventdata, handles,trialData,varargin)
% hObject    handle to FitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

currentModel = getappdata(parent,'CurrentModel');
currentDataSet = getappdata(parent,'CurrentDataSet');

selected = unique(cat(2,get(handles.modelList,'value')-1,currentModel));

selected(selected==0) = [];

modelData = getappdata(parent,'modelData');
regData = getappdata(parent,'regData');
if nargin < 4 || isempty(trialData)
    trialData = getappdata(parent,'trialData');
    trialData = trialData(currentDataSet);
end

use_parallel = true; %Uses parallel computing toolbox if available (not yet implemented)

for i = 1:length(selected)
    currmodel = selected(i);
    regs = ismember([regData(currentDataSet).regressors.code],modelData(currentDataSet).models(currmodel).regressors);
    
    R = regData(currentDataSet).regressors(regs);   
    binvolume = regData(currentDataSet).regressors(3).value;
%     dlg = msgbox('Fitting...');
      fprintf('\nFitting model %s...',modelData(currentDataSet).models(currmodel).label);
    fit = modelFit(trialData,R, 'fullonly',~get(handles.subModelCheck,'value'),...
                    'Firth',get(handles.firthCheck,'value'),...
                    'multassign',get(handles.multiBinCheck,'value'),...
                    'regularization',modelData(currentDataSet).models(currmodel).Hreg,...
                    'binvolume',binvolume,varargin{:});
%     if ishandle(dlg)
%         delete(dlg)            
%     end
    modelData(currentDataSet).models(selected(i)).fit = fit;
    setappdata(parent,'modelData',modelData);
    
end
   


    


%--------------------------------

function Update(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');

currentModel = getappdata(parent,'CurrentModel');
% currentRegressorGroup = getappdata(parent,'currentRegressorGroup');
currentDataSet = getappdata(parent,'CurrentDataSet');

modelData = getappdata(parent,'modelData');
regData = getappdata(parent,'regData');

if currentDataSet == 0
    fprintf('\nNo Data Set selected...')
    return
end


if currentDataSet > length(modelData)
   modelData(currentDataSet) = makeModelData([]);
    modelData(currentDataSet).models(:) = [];
end

if currentModel > length(modelData(currentDataSet).models);
    currentModel = length(modelData(currentDataSet).models);
    setappdata(parent,'CurrentModel',currentModel);
    set(handles.modelList,'value',currentModel+1)

end

if currentModel > 0
    modelData(currentDataSet).models(currentModel).regressors = unique(modelData(currentDataSet).models(currentModel).regressors);
    regcodes = [regData(currentDataSet).regressors.code];
    [q,rgindx] = ismember(modelData(currentDataSet).models(currentModel).regressors,regcodes);
    if any(rgindx)
        modelData(currentDataSet).models(currentModel).regLabels = {regData(currentDataSet).regressors(rgindx).label};
        setappdata(parent,'modelData',modelData);
    end
end

if isnumeric(modelData(currentDataSet).models)
    modelData(currentDataSet) = makeModelData;
end

liststr  = {'New Model'};
set(handles.modelList,'string',cat(2,liststr,{modelData(currentDataSet).models.label}));

liststr  = {'Add Selected RegressorGroup Groups'};
if currentModel > length(modelData(currentDataSet).models) || isempty(currentModel ) || currentModel == 0 
    currentModel = 0;
    setappdata(parent,'CurrentModel',currentModel);
    set(handles.modelList,'value',1);
end

if ~isempty(currentModel) && currentModel > 0 
    regcodes = [regData(currentDataSet).regressors.code];
    [q,rgindx] = ismember(modelData(currentDataSet).models(currentModel).regressors,regcodes);

%     modelregs =modelData(currentDataSet).models(currentModel).regressors;
%     a = ismember([regData(currentDataSet).regressors.code],modelregs);
    lbl = strcat(cellfun(@num2str,num2cell(1:sum(q)),'uniformoutput',0),{'. '},{regData(currentDataSet).regressors(rgindx).label});
    val = get(handles.regressorGroupList,'value');
    if val > length(lbl)+1
         set(handles.regressorGroupList,'value',length(lbl)+1);
    end    
    
    set(handles.regressorGroupList,'string',cat(2,liststr,lbl));
%     set(handles.regressorGroupList,'string',cat(2,liststr,{regData(currentDataSet).regressors(ismember([regData(currentDataSet).regressors.code],modelData(currentDataSet).models(currentModel).regressors)).label}));
    set(handles.label,'string',modelData(currentDataSet).models(currentModel).label);
else
    set(handles.label,'string','none selected');
end



% --- Executes on button press in subModelCheck.
function subModelCheck_Callback(hObject, eventdata, handles)
% hObject    handle to subModelCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subModelCheck


% % --- Executes during object creation, after setting all properties.
% function regressorGroupList_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to regressorGroupList (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: listbox controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in multiBinCheck.
function multiBinCheck_Callback(hObject, eventdata, handles)
% hObject    handle to multiBinCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiBinCheck


% --- Executes on button press in firthCheck.
function firthCheck_Callback(hObject, eventdata, handles)
% hObject    handle to firthCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of firthCheck


% --------------------------------------------------------------------
function ModelListContext_Callback(hObject, eventdata, handles)
% hObject    handle to ModelListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function addRegressorsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to addRegressorsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
rm = getappdata(parent,'RegManager');
currentDataSet= getappdata(parent,'CurrentDataSet');
currentModel= getappdata(parent,'CurrentModel');

if isempty(rm) || ~ishandle(rm) || isequal(currentModel,0) || isequal(currentDataSet,0)
    return
end

rmh = guidata(rm);

regData = getappdata(parent,'regData');
modelData = getappdata(parent,'modelData');

selectedrg = get(rmh.regressorGroupList,'value') - 1;
selectedrg(selectedrg==0) = [];

modelData(currentDataSet).models(currentModel).regressors =...
    sort(cat(2,modelData(currentDataSet).models(currentModel).regressors,regData(currentDataSet).regressors(selectedrg).code));

setappdata(parent,'modelData',modelData);

                                                            
 Update(hObject, eventdata, handles)

% --------------------------------------------------------------------
function duplicateMenu_Callback(hObject, eventdata, handles)
% hObject    handle to duplicateMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');
modelData = getappdata(parent,'modelData');
currentDataSet= getappdata(parent,'CurrentDataSet');
currentModel= getappdata(parent,'CurrentModel');


modelData(currentDataSet).models(end+1)= modelData(currentDataSet).models(currentModel);
modelData(currentDataSet).models(end).code = modelData(currentDataSet).codeincr + 1;
modelData(currentDataSet).models(end).label =  [modelData(currentDataSet).models(end).label,' cpy'];
modelData(currentDataSet).codeincr = modelData(currentDataSet).codeincr + 1;


setappdata(parent,'modelData',modelData);

Update(hObject, eventdata, handles)


