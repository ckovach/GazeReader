function varargout = TrialManager(varargin)
% TRIALMANAGER M-file for TrialManager.fig
%      TRIALMANAGER, by itself, creates a new TRIALMANAGER or raises the existing
%      singleton*.
%
%      H = TRIALMANAGER returns the handle to a new TRIALMANAGER or the handle to
%      the existing singleton*.
%
%      TRIALMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIALMANAGER.M with the given input arguments.
%
%      TRIALMANAGER('Property','Value',...) creates a new TRIALMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrialManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrialManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrialManager

% Last Modified by GUIDE v2.5 01-Jan-2008 15:07:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrialManager_OpeningFcn, ...
                   'gui_OutputFcn',  @TrialManager_OutputFcn, ...
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


% --- Executes just before TrialManager is made visible.
function TrialManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrialManager (see VARARGIN)

% Choose default command line output for TrialManager
handles.output = hObject;

parent = varargin{1};
setappdata(handles.figure1,'parent',parent );

ch = getappdata(parent,'children');
ch(end+1) = handles.figure1;
setappdata(parent,'children',ch);

% Update handles structure
guidata(hObject, handles);

% trialData = getappdata(parent,'trialData');

% if isempty(trialData)
%     setapppdata(parent,'trialData',makeTrialData([]));
% end

tmfuns.update = @(varargin) Update([],[],handles,varargin{:});
tmfuns.updateFixations = @(varargin)UpdateFixations([], [], handles,varargin{:});
tmfuns.updateAllDataSets =  @(varargin)UpdateAllDataSets([], [], handles,varargin{:});

tmfuns.get_time_intvl_data = @() get_time_intvl_data(hObject,[],handles);

setappdata(parent,'trialManagerFunctions',tmfuns);

crdat = getappdata(parent,'CurrentDataSet');
etd = getappdata(parent,'eyetrackerHeaderData');
% 
% for i = 1:length(etd)
%     setappdata(parent,'CurrentDataSet',i)
%     Update(hObject, eventdata, handles) 
%     UpdateFixations(hObject, eventdata, handles)
% end

% for i = 1:length(etd)
    setappdata(parent,'CurrentDataSet',crdat)
    Update(hObject, eventdata, handles) 
    UpdateFixations(hObject, eventdata, handles)
% end

setappdata(parent,'CurrentDataSet',crdat)
 

% UIWAIT makes TrialManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% parent = getappdata(handles.figure1,'parent');
% 
% currentTrial = getappdata(parent,'CurrentTrial');
% trialData = getappdata(parent,'trialData');



% --- Outputs from this function are returned to the command line.
function varargout = TrialManager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function TrialNumber_Callback(hObject, eventdata, handles)
% hObject    handle to TrialNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrialNumber as text
%        str2double(get(hObject,'String')) returns contents of TrialNumber as a double

trnum = str2double(get(handles.TrialNumber,'string'));

if ~isempty(trnum)
    parent = getappdata(handles.figure1,'parent');

    setappdata(parent,'CurrentTrial',trnum);
    set(handles.trialList,'value',trnum+1) %,'value',currentTrial+1)
end

Update(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function TrialNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrialNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trialList.
function trialList_Callback(hObject, eventdata, handles)
% hObject    handle to trialList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns trialList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trialList


parent = getappdata(handles.figure1,'parent');

selectedTrial = get(handles.trialList,'value')-1;

setappdata(parent,'CurrentTrial',selectedTrial(1));
setappdata(parent,'CurrentBin',0);
setappdata(parent,'CurrentFixation',0);

Update(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function trialList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function buildFromDataSet_Callback(hObject, eventdata, handles)
% hObject    handle to buildFromDataSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

ch = getappdata(parent,'children');
% ch(end+1) = TrialBuilder(parent);
ch(end+1) = EventManager(parent);
setappdata(parent,'children',ch);


% eventData = getappdata(parent,'eventData');
% 
% evtcodenum = cellfun(@num2str,{eventData.xdat.id},'uniformoutput',0);
% evtcodes = cellfun(@num2str,eventData.codes([eventData.xdat.id]),'uniformoutput',0);
% cstrlen = cellfun(@length,evtcodenum);
% spaces1 = cellfun(@(n) repmat(' ',1,n),mat2cell(10-2*cstrlen,1,ones(1,length(cstrlen))),'UniformOutput',false);
% evt = [eventData.xdat.startT];
% evt = evt-evt(1);
% evttimes = cellfun(@num2str,mat2cell(evt,1,ones(1,length(evt))),'uniformoutput',0);
% tstrlen = cellfun(@length,evttimes);
% spaces2 = cellfun(@(n) repmat(' ',1,n),mat2cell(15-2*tstrlen,1,ones(1,length(cstrlen))),'UniformOutput',false);
% liststr = strcat(evttimes,spaces2,evtcodenum,spaces1,evtcodes);
% 
% fg = figure;
% set(fg,'MenuBar','none');
% lb = uicontrol(fg,'String',liststr,'Style','Listbox','units','normalized','string',liststr,'position',[.1 .1 .8 .8]);

% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function assignSelectedImagesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to assignSelectedImagesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
imh = getappdata(parent,'ImageManager');
if ~ishandle(imh)
    return
end
imhandles = guidata(imh);
trialData = getappdata(parent,'trialData');
imageData = getappdata(parent,'imageData');

activeIM = get(imhandles.imageList,'Value')-1;
activeIM(activeIM== 0) = [];
activeTrials = get(handles.trialList,'Value')-1;
activeTrials (activeTrials == 0) = [];

currentDataSet = getappdata(parent,'CurrentDataSet');
% m2c = @(x) mat2cell(x,1,ones(size(x)));

if isempty(activeIM) 
    empty = repmat({[]},1,length(activeTrials));
    [trialData(currentDataSet).trials(activeTrials).image] = empty{:};
else
    if length(activeIM) == 1
        activeIM = repmat(activeIM,1,length(activeTrials));

    elseif length(activeIM) ~= length(activeTrials)
        errordlg('Number of selected images and selected trials mismatch. Try again.')
        return
    end

    [trialData(currentDataSet).trials(activeTrials).image] = imageData.images(activeIM).code;
    [imageData.images(activeIM).activeTrials] = trialData(currentDataSet).trials(activeTrials).number;
    setappdata(parent,'imageData',imageData)
end

setappdata(parent,'trialData',trialData)
Update(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function assignSelectedBinsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to assignSelectedBinsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parent = getappdata(handles.figure1,'parent');

bmhandles = guidata(getappdata(parent,'BinManager'));
trialData = getappdata(parent,'trialData');
binData = getappdata(parent,'binData');
currentDataSet = getappdata(parent,'CurrentDataSet');

activeBG = get(bmhandles.binGroupList,'Value')-1;
activeBG(activeBG == 0) = [];
activeTrials = get(handles.trialList,'Value')-1;
activeTrials (activeTrials == 0) = [];


% m2c = @(x) mat2cell(x,1,ones(size(x)));

% if length(activeBG) == 1
%     activeBG = repmat(activeBG,1,length(activeTrials));
% elseif length(activeBG) ~= length(activeTrials)
%     error('Number of selected bin groups and selected trials mismatch. Try again.')
% end

bins = repmat({[binData.groups(activeBG).code]},1,length(activeTrials));
trials = repmat({[trialData(currentDataSet).trials(activeTrials).code]},1,length(activeBG));

% [trialData(currentDataSet).trials(activeTrials).binGroup] = binData.groups(activeBG).code;
% [binData.groups(activeBG).activeTrials] = trialData(currentDataSet).trials(activeTrials).number;
[trialData(currentDataSet).trials(activeTrials).binGroup] = bins{:};
[binData.groups(activeBG).activeTrials] = trials{:};

setappdata(parent,'trialData',trialData)
setappdata(parent,'binData',binData)

Update(hObject, eventdata, handles)
UpdateFixations(hObject, eventdata, handles)


% --------------------------------------------------------------------
function assignDuplicateBinsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to assignDuplicateBinsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

bmhandles = guidata(getappdata(parent,'BinManager'));
trialData = getappdata(parent,'trialData');
binData = getappdata(parent,'binData');
currentDataSet = getappdata(parent,'CurrentDataSet');

activeBG = get(bmhandles.binGroupList,'Value')-1;
activeBG(activeBG == 0) = [];
activeTrials = get(handles.trialList,'Value')-1;
activeTrials (activeTrials == 0) = [];


newBinData = binData;
newBinData.groups = repmat(binData.groups(activeBG),1,length(activeTrials));
codec = mat2cell(1:length(newBinData.groups),1,ones(1,length(newBinData.groups)));
[newBinData.groups.code] = codec{:};
% bins = repmat({[binData.groups(activeBG).code]},1,length(activeTrials));
% trials = repmat({[trialData(currentDataSet).trials(activeTrials).number]},1,length(activeBG));

newbinindex = length(binData.groups) + repmat((1:length(activeBG))-1,1,length(activeTrials)) + kron(1:length(activeTrials),ones(1,length(activeBG)));

binData = makeBinData(binData,newBinData);
newcodes = mat2cell([binData.groups(newbinindex).code],1,length(activeBG)*ones(1,length(activeTrials)));

[trialData(currentDataSet).trials(activeTrials).binGroup] = newcodes{:};

acttr = kron([trialData(currentDataSet).trials(activeTrials).number],ones(1,length(activeBG)));
acttrc = mat2cell(acttr,1,ones(1,length(acttr)));
newlabels = strcat({binData.groups(newbinindex).label},{'_Trial_'},cellfun(@num2str,acttrc,'UniformOutput',0));

[binData.groups(newbinindex).activeTrials] = acttrc{:};
[binData.groups(newbinindex).label] = newlabels{:};


setappdata(parent,'trialData',trialData)
setappdata(parent,'binData',binData)

bmfun = getappdata(parent,'binManagerFunctions');

Update(hObject, eventdata, handles)
UpdateFixations(hObject, eventdata, handles)

bmfun.update(parent);


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Update(hObject, eventdata, handles, varargin)

if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
else
    parent = varargin{1};
end

currentTrial = getappdata(parent,'CurrentTrial');
trialData = getappdata(parent,'trialData');
% expEventData = getappdata(parent,'expEventData');
imageData = getappdata(parent,'imageData');
binData = getappdata(parent,'binData');
currentDataSet = getappdata(parent,'CurrentDataSet');

if currentDataSet == 0
    fprintf('No data set selected...\n')
    return
end
if isempty(currentTrial)
    currentTrial = 0;
end
% if ishandle(handles.figure1)
    if ~isempty([trialData(currentDataSet).trials.code])
        
%         emptycell
        trnum = cellfun(@num2str,{trialData(currentDataSet).trials.number},'uniformoutput',0);
        trlen = cellfun(@length,trnum);
        spaces1 = cellfun(@(n) repmat(' ',1,n),mat2cell(8-2*trlen,1,ones(1,length(trlen))),'UniformOutput',false);
        
        startTime = cellfun(@num2str,{trialData(currentDataSet).trials.startTime},'uniformoutput',0);
        sttlen = cellfun(@length,startTime);
        spaces2 = cellfun(@(n) repmat(' ',1,n),mat2cell(20-2*sttlen,1,ones(1,length(sttlen))),'UniformOutput',false);

        dur = [trialData(currentDataSet).trials.stopTime] - [trialData(currentDataSet).trials.startTime];
        if ~isempty(dur)
        duration = cellfun(@num2str,mat2cell(dur,1,ones(1,length(dur))),'uniformoutput',0);
        durlen = cellfun(@length,duration);
        spaces4 = cellfun(@(n) repmat(' ',1,n),mat2cell(20-2*durlen ,1,ones(1,length(durlen))),'UniformOutput',false);
        % stopTime = cellfun(@num2str,{trialData(currentDataSet).trials.stopTime},'uniformoutput',0);
        % stoptlen = cellfun(@length,stopTime );
        % spaces4 = cellfun(@(n) repmat(' ',1,n),mat2cell(20-2*stoptlen ,1,ones(1,length(stoptlen))),'UniformOutput',false);

    %     startcodes = cellfun(@num2str,expEventData.codes([trialData(currentDataSet).trials.startCode]),'uniformoutput',0);
        startcodes = {trialData(currentDataSet).trials.startCodeLabel};
        startlen = cellfun(@length,startcodes);
        spaces3 = cellfun(@(n) repmat(' ',1,n),mat2cell(30-2*startlen,1,ones(1,length(startlen))),'UniformOutput',false);

    %     stopcodes = cellfun(@num2str,expEventData.codes([trialData(currentDataSet).trials.stopCode]),'uniformoutput',0);
        stopcodes = {trialData(currentDataSet).trials.stopCodeLabel};
        stoplen = cellfun(@length,stopcodes);
        spaces5 = cellfun(@(n) repmat(' ',1,n),mat2cell(30-2*stoplen,1,ones(1,length(stoplen))),'UniformOutput',false);


        trliststr = strcat(trnum,spaces1,startTime,spaces2,duration,spaces4,startcodes,spaces3,stopcodes,spaces5);
    else
        trliststr = {};
    end

    trliststr = [{'#     StartT       Duration       StartCode             StopCode    '},trliststr ];

    set(handles.trialList,'string',trliststr) %,'value',currentTrial+1)
    if max(get(handles.trialList,'value')) > length(trliststr)
       set(handles.trialList,'value',1) 
    end
    
    while currentTrial > length([trialData(currentDataSet).trials])
        currentTrial  = currentTrial -1;
    end
    setappdata(parent,'CurrentTrial',currentTrial);

    if currentTrial ~=0 &&  ~isempty([trialData(currentDataSet).trials.number])

        currtr = trialData(currentDataSet).trials(currentTrial);

        set(handles.startTimeText,'string',num2str(currtr.startTime));
        set(handles.stopTimeText,'string',num2str(currtr.stopTime));   
        set(handles.startCodeStringText,'string',currtr.startCodeLabel);
    %     set(handles.startCodeStringText,'string',expEventData.codes(currtr.startCode));
        set(handles.stopCodeStringText,'string',currtr.stopCodeLabel);
    %     set(handles.stopCodeStringText,'string',expEventData.codes(currtr.stopCode));
    
        if isnumeric(currtr.startCode), str = num2str(currtr.startCode); else  str = currtr.startCode;end
        set(handles.startCodeText,'string',str);
        if isnumeric(currtr.startCode), str = num2str(currtr.stopCode); else  str = currtr.stopCode;end
        set(handles.stopCodeText,'string',str);
        if ~isempty(currtr.image)
            trim = ismember([imageData.images.code],currtr.image);
            set(handles.activeImageText,'string',imageData.images(trim).filename);
        else
            set(handles.activeImageText,'string','none');
        end
        if ~isempty(currtr.binGroup)   
            trbn = ismember([binData.groups.code],currtr.binGroup);
            set(handles.binList,'string',{binData.groups(trbn).label} );
        else
            set(handles.binList,'string','none active');
            set(handles.binList,'value',1);
        end

    else
        set(handles.startTimeText,'string','');
        set(handles.stopTimeText,'string','');
        set(handles.startCodeStringText,'string','');
        set(handles.stopCodeStringText,'string','');
        set(handles.startCodeText,'string','');
        set(handles.stopCodeText,'string','');

        set(handles.activeImageText,'string','none');
        set(handles.binList,'string','none');
    end
end

if currentTrial ~=0
    
    imh = getappdata(parent,'ImageManager');
    bmh = getappdata(parent,'BinManager');
    dmh = getappdata(parent,'DataDisplayManager');
    rmh = getappdata(parent,'RegManager');
    
%     if ~isempty(imh)
    if ishandle(imh)
        imcodes = [imageData.images.code];
        
        phandles = guidata(parent);
        imfun= getappdata(parent,'ImageManagerFunctions');
        activeim =  find(ismember(imcodes,trialData(currentDataSet).trials(currentTrial).image));
        if ~isempty(activeim ) 
            if  ishandle(imh)
                imhandles = guidata(imh);        
                set(imhandles.imageList,'value',activeim+1);
            end
            setappdata(parent,'CurrentImage',activeim(1));
            imfun.UpdateImage(parent);
        else
            cla(phandles.axes1)
        end
    end
    
%     if ~isempty(bmh)
    if ishandle(bmh)
        bmcodes = [binData.groups.code];
        bmfun = getappdata(parent,'binManagerFunctions');
        activebn =  find(ismember(bmcodes,trialData(currentDataSet).trials(currentTrial).binGroup));
        if ~isempty(activebn)
%             if ishandle(bmh)
                bmhandles = guidata(bmh);
                set(bmhandles.binGroupList,'value',activebn+1);
%             end
            setappdata(parent,'CurrentBinGroup',activebn(1));
%             bmfun.update();
            if get(bmhandles.displayCheck,'value')
                bmfun.clear(parent);
                bmfun.draw(parent);
            end
        else
            bmfun.clear(parent);
        end
    end    
    
	if ~isempty(dmh) && ishandle(dmh) && ~isempty(hObject)
        try
            dmfun = getappdata(parent,'DataDisplayFunctions');
            dmfun.update();
        catch
            warning('Failed to update display')
        end
    end
    if ~isempty(rmh) && ishandle(rmh)
        try
            rmfun = getappdata(parent,'RegManagerFunctions');
            rmfun.update();
        catch
        end
    end
end

if ishandle(handles.figure1)
    blstr = get(handles.binList,'string');
    blv = get(handles.binList,'value');
    if blv >  length(blstr)
        set(handles.binList,'value',1);
    end



    set(handles.TrialNumber,'string',num2str(currentTrial)) %,'value',currentTrial+1)
    if currentTrial~=0
        set(handles.idCodeEdit,'string',num2str(trialData(currentDataSet).trials(currentTrial).code)) %,'value',currentTrial+1)
    else
        set(handles.idCodeEdit,'string','0') %,'value',currentTrial+1)
    end    

    if ~isempty(hObject)
        activefigure(handles.figure1)
    end
end
if  currentTrial ~= 0 && ~isempty([trialData(currentDataSet).trials(currentTrial).fixations]) 
    setappdata(parent,'CurrentFixation',trialData(currentDataSet).trials(currentTrial).fixations(1))
else
    setappdata(parent,'CurrentFixation',[])
end

% --------------------------------------------------------------------
function removeBinMenu_Callback(hObject, eventdata, handles)
% hObject    handle to removeBinMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
trialData = getappdata(parent,'trialData');
currentTrial = getappdata(parent,'CurrentTrial');
currentDataSet = getappdata(parent,'CurrentDataSet');

activeTrials = get(handles.trialList,'value')-1;
activeTrials(activeTrials==0) = [];

if isempty(trialData(currentDataSet).trials(currentTrial).binGroup);
    return
end
selectedBins = get(handles.binList,'value');
activeBins = trialData(currentDataSet).trials(currentTrial).binGroup(selectedBins);


for i = 1:length(activeTrials)
    bg = trialData(currentDataSet).trials(activeTrials(i)).binGroup;    
    trialData(currentDataSet).trials(activeTrials(i)).binGroup(ismember(bg,activeBins)) = [];
end

% trialData(currentDataSet).trials(currentTrial).binGroup(activeBin) = [];

setappdata(parent,'trialData',trialData);

if selectedBins(1) > 1
    set(handles.binList,'value',selectedBins(1) - 1);
end
   
Update(hObject, eventdata, handles);
UpdateFixations(hObject, eventdata, handles)

% --------------------------------------------------------------------
function BinListContextMenu_Callback(hObject, eventdata, handles)
% hObject    handle to BinListContextMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function trialListContextMenu_Callback(hObject, eventdata, handles)
% hObject    handle to trialListContextMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function removeTrialMenu_Callback(hObject, eventdata, handles)
% hObject    handle to removeTrialMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

trialData = getappdata(parent,'trialData');
expEventData = getappdata(parent,'expEventData');
currentDataSet = getappdata(parent,'CurrentDataSet');

selected = get(handles.trialList,'value')-1;

selected(selected==0) = [];


if ~isempty(expEventData)
    resetevts = ismember([expEventData(currentDataSet).events.code],[trialData(currentDataSet).trials(selected).startEventCode,...
                trialData(currentDataSet).trials(selected).stopEventCode]);
    c = mat2cell(zeros(1,sum(resetevts)),1,ones(1,sum(resetevts)));
    [expEventData(currentDataSet).events(resetevts).type] = c{:};
end

trialData(currentDataSet).trials(selected) = [];

% expEventData.trialOnsets(selected) = [];
% expEventData.trialOnsetCodes(selected) = [];
% expEventData.trialEnds(selected) = [];
% expEventData.trialEndCodes(selected) = [];

if ~isempty(trialData(currentDataSet).trials)
    
    trnum = mat2cell(1:length(trialData(currentDataSet).trials),1,ones(1,length(trialData(currentDataSet).trials)));
    [trialData(currentDataSet).trials.number] = trnum{:};

end

setappdata(parent,'trialData',trialData);
setappdata(parent,'expEventData',expEventData);
setappdata(parent,'CurrentTrial',selected(1)-1);
set(handles.trialList,'value',selected(1));
Update(hObject, eventdata, handles);


%------------------------------------------------
function idCodeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to idCodeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of idCodeEdit as text
%        str2double(get(hObject,'String')) returns contents of idCodeEdit as a double

trcode = str2double(get(handles.TrialNumber,'string'));
parent = getappdata(handles.figure1,'parent');
currentDataSet = getappdata(parent,'CurrentDataSet');

trialData = getappdata(parent,'trialData');
trcodes = [trialData(currentDataSet).trials.code];

trnum = find(trcodes == trcode);
 
if ~isempty(trnum)

    setappdata(parent,'CurrentTrial',trnum);
    set(handles.trialList,'value',trnum+1) %,'value',currentTrial+1)
end

Update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function idCodeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to idCodeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%------------------------------------------------
function UpdateFixations(hObject, eventdata, handles,varargin)

%Recomputes fixation bin membership based on current information

fixation_position_field = 'meanPos';
default_etfs  = 1; % Default sampling rate if no other is specified

if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
else
    parent = varargin{1};
end
  
trialData = getappdata(parent,'trialData');

currentDataSet = getappdata(parent,'CurrentDataSet');
binData = getappdata(parent,'binData');

if isempty([trialData(currentDataSet).trials.number]) || isempty(binData)
    return
end

fixdat = getappdata(parent,'fixationData');
% headerData = getappdata(parent,'eyetrackerHeaderData');
% expEventData = getappdata(parent,'eyetrackerHeaderData');

% headerData =headerData(currentDataSet);
if isempty(fixdat)
    return
end

fixdat = fixdat(currentDataSet);
fxstarts = [fixdat.fix.startT];

trdat = trialData(currentDataSet);

% currentTrial = getappdata(parent,'CurrentTrial');
% selected = getappdata(handles.trialManager,'value')-1;
% selected(selected==0) = [];
% updateTrials = unique([currentTrial,selected]);


% if isfield(expEventData ,'xdat') && ~isempty(expEventData.xdat)
%     sttime = expEventData.xdat(1).startT;
% else
%     sttime = 0;
% end

CurrentDataSet = getappdata(parent,'CurrentDataSet');

bingroupcodes = [binData.groups.code];


screenData = getappdata(parent,'screenData');
if CurrentDataSet > length(screenData)
    screenData(CurrentDataSet) = screenData(1);
end

screenres = screenData(CurrentDataSet).res;
if isfield(fixdat,'etfs') && ~isempty(fixdat.etfs) && ~isequal(fixdat.etfs,0)
    etfs = fixdat.etfs;

else
    etfs = default_etfs;
end

for i = 1:length(trdat.trials)  %Eventually change to modify only altered trials
    
    trdat.trials(i).samplePts = floor([trdat.trials(i).startTime:1./etfs:trdat.trials(i).stopTime]*etfs);
    trdat.trials(i).seg = currentDataSet;
    if isempty(trdat.trials(i).startTime) || isempty(trdat.trials(i).stopTime)
        trdat.trials(i).startTime=Inf;
        trdat.trials(i).stopTime = Inf;
    end
    if isempty(trdat.trials(i).code)
        trdat.trials(i).code = max([trdat.trials.code])+1;
        trdat.codeincr = trdat.trials(i).code;
    end
    if isempty(trdat.trials(i).number)
        trdat.trials(i).number = nan;
    end  
    fixInTrial = find( fxstarts >= trdat.trials(i).startTime & fxstarts <= trdat.trials(i).stopTime);
    trdat.trials(i).fixations = fixInTrial;
    trdat.trials(i).fixOnsetTimes = fxstarts( fixInTrial ); 
    trdat.trials(i).fixTrialTimes = fxstarts( fixInTrial ) - trdat.trials(i).startTime;
    
    trialbingrs = binData.groups(ismember(bingroupcodes, trdat.trials(i).binGroup));
    
    nbins = cumsum([trialbingrs.nbin]);
    nfix = length(fixInTrial);
    
    if isempty(nbins)
        nbins = 0;
    end
    
    fixlscr = cat(1,fixdat.fix(fixInTrial).(fixation_position_field));
    
    if nfix > 0
        fixlocs = repmat({{fixlscr./repmat(screenres,nfix,1)}},1,length(trialbingrs));
    end
    
    binmemberfuns = {trialbingrs.isinside}; %functions for deciding bin memberhsip
    
    if ~isempty(trialbingrs)
        ccbingr = mat2cell(mat2cell(trialbingrs,1,ones(1,length(trialbingrs))),1,ones(1,length(trialbingrs)));
        multfun = @(fun,arg1, arg2) cellfun(fun,arg1,arg2,'uniformoutput',0);
        if nfix > 0
            binmembership = cellfun(multfun, binmemberfuns, ccbingr , fixlocs); 
        else
            binmembership = {spalloc(0,nbins(end),0)};            
        end
        
        binareafuns = {trialbingrs.area}; %functions that return bin area
        multfun2 = @(fun,arg) cellfun(fun,arg,'uniformoutput',0);
          binareas = cellfun(multfun2, binareafuns, ccbingr ); 
          
         trdat.trials(i).binareas = cat(1,binareas{:});
%         for bg = 1:length(trialbingrs)
%             binmembership{bg} = binmemberfuns{gb}(ccbingr{bg}{1}, fixlocs{bg}{1});
%         end
        binmembership = cat(2,binmembership{:});
    else    
        binmembership =sparse(zeros(nfix,0));
        trdat.trials(i).binareas = [];       
    end
    
    trdat.trials(i).nfix = nfix;
    trdat.trials(i).nbin = nbins(end);
    [ir,jc] = sparseij(binmembership');
    
    djc = diff(jc);
    if isempty(djc)
        djc = 0;
    end
    
%      djc(jc(1:end-1) == jc(end)) = [];
     jc(jc == jc(end)) = [];
%     trdat.trials(i).fixbin = ir(jc(1:end-1)+1)+1;
    trdat.trials(i).fixmat = binmembership;
    if nfix == 0
        continue
    end
    fxb =  ir(jc(1:end)+1)+1;
    fxb(djc==0) = 0;
    
    trdat.trials(i).fixbin = fxb(:);
    
    
    if any(djc  > 1)  %Check for multiple assignment
         
        switch binData.precedence

            case {'group_order','array_order'}  % Simply assign to the first bin, which is already done
            case 'nearest_to_center'
      
                bincenters = cat(1,trialbingrs.centers);
                
                for fx = find(djc>1)'
                    
                    fixloc = fixlscr(fx,:);
                    fulljc = [jc; length(ir)];
                    ovlapbins = ir( (fulljc(fx):fulljc(fx+1)-1)+1)+1; %bins in multiple assignment
                    
                    d = sqrt(sum((repmat(fixloc,size(ovlapbins,1),1) - bincenters(ovlapbins,:).*repmat(screenres,length(ovlapbins),1)).^2,2));
                    [md,mindex] = min(d);
                    
                    trdat.trials(i).fixbin(fx) = ovlapbins(mindex);
%                     trdat.trials(i).fixbin(jc==0) = 0;
                    
                end
            otherwise
                error('Unrecognized multiple assignment code')               

        end
    end
end

trialData(currentDataSet) = trdat;

setappdata(parent,'trialData',trialData)

%------- Updates fixations for all data sets
function UpdateAllDataSets(hObject, eventdata, handles,varargin)
    

if ishandle(handles.figure1)
    parent = getappdata(handles.figure1,'parent');
else
    parent = varargin{1};
end

trialData = getappdata(parent,'trialData');
 
for i = 1:length(trialData)
    setappdata(parent,'CurrentDataSet',i)
     UpdateFixations(hObject, eventdata, handles,varargin)
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
if ishandle(parent)
    setappdata(parent,'activeControl','main');
end

% function to get time intervals 
function intvl_data = get_time_intvl_data( hObject, eventdata,handles)

parent = getappdata(handles.figure1,'parent');

trialData = getappdata(parent, 'trialData');
fixationData = getappdata(parent, 'fixationData');
CurrentDataSet = getappdata(parent, 'CurrentDataSet');
time_intvl = getappdata(parent,'TSampIntvl');


for i = 1:length(trialData(CurrentDataSet).trials)
    fxindxs = trialData(CurrentDataSet).trials(i).fixations;
	if    trialData(CurrentDataSet).trials(i).nfix == 0
        fxindxs = max([trialData(CurrentDataSet).trials(1:i).fixations]);
    elseif  fxindxs(1) == 1;
        fxindxs = [fxindxs(1),fxindxs];
    else
        fxindxs = [fxindxs(1)-1,fxindxs];
    end        
    fxonset = double([fixationData(CurrentDataSet).fix(fxindxs).startT]);
    fxbins = ceil((fxonset(2:end)-trialData(CurrentDataSet).trials(i).startTime)./time_intvl);
    if ~isempty(fxbins)
        fxnbins = [fxbins(1),diff(fxbins)]';
    else
        fxnbins = [];
    end
    
    intvl_data(i).times = double((trialData(CurrentDataSet).trials(i).startTime:time_intvl:trialData(CurrentDataSet).trials(i).stopTime)' + time_intvl./2);            
    intvlt =  intvl_data(i).times -trialData(CurrentDataSet).trials(i).startTime;
    fxi = [0;zeros(size(intvl_data(i).times))];
    fxi(fxbins+1) = 1:trialData(CurrentDataSet).trials(i).nfix;   
    intvl_data(i).fxi = fxi(2:end);
    intvl_data(i).trtimes = intvlt;
    intvl_data(i).fxnbins = fxnbins;
%     fxt =diff(ceil(fixrt./time_intvl));
    
    fxi(1) = 1;
    q = zeros(size(fxi));
    q(fxi~=0) = fxonset(1:end)-[0,fxonset(1:end-1)];
    q = cumsum(q);
    
    intvl_data(i).fxtimes = intvl_data(i).times-q(2:end);
end



