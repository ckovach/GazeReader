function varargout = RegManager(varargin)
% REGMANAGER M-file for RegManager.fig
%      REGMANAGER, by itself, creates a new REGMANAGER or raises the existing
%      singleton*.
%
%      H = REGMANAGER returns the handle to a new REGMANAGER or the handle to
%      the existing singleton*.
%
%      REGMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGMANAGER.M with the given input arguments.
%
%      REGMANAGER('Property','Value',...) creates a new REGMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RegManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RegManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RegManager

% Last Modified by GUIDE v2.5 06-Feb-2008 19:41:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RegManager_OpeningFcn, ...
                   'gui_OutputFcn',  @RegManager_OutputFcn, ...
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


% --- Executes just before RegManager is made visible.
function RegManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RegManager (see VARARGIN)


parent = varargin{1};

setappdata(handles.figure1,'parent',parent);
setappdata(parent,'RegManager',handles.figure1);

children = getappdata(parent,'children');
children(end+1) = handles.figure1;
setappdata(parent,'children',children)

setappdata(parent,'CurrentRegressorGroup',0);
setappdata(parent,'CurrentRegressor',0);

setappdata(parent,'CurrentModel',1);

rmfuns.update = @() UpdateFields([],[],handles);
rmfuns.designMatrix= @(varargin) designMatrix([],[],handles,varargin{:});


rmfuns.reset= @() resetRegressors(hObject,[],handles);

setappdata(parent,'RegManagerFunctions',rmfuns);

InitializeRegressors(hObject, eventdata, handles)

UpdateFields(hObject, eventdata, handles)

% Choose default command line output for RegManager
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RegManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RegManager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in regressorGroupList.
function regressorGroupList_Callback(hObject, eventdata, handles)
% hObject    handle to regressorGroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns regressorGroupList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from regressorGroupList


parent = getappdata(handles.figure1,'parent');

selected = get(handles.regressorGroupList,'value');

setappdata(parent,'CurrentRegressorGroup',selected(1) - 1)

UpdateFields(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function regressorGroupList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regressorGroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in regList.
function regList_Callback(hObject, eventdata, handles)
% hObject    handle to regList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns regList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from regList

parent = getappdata(handles.figure1,'parent');

selected = get(handles.regList,'value');

setappdata(parent,'CurrentRegressor',selected - 1)

UpdateFields(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function regList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function importFromWorkspaceMenu_Callback(hObject, eventdata, handles)
% hObject    handle to importFromWorkspaceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function importFromTextMenu_Callback(hObject, eventdata, handles)
% hObject    handle to importFromTextMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function regGroupRemoveMenu_Callback(hObject, eventdata, handles)
% hObject    handle to regGroupRemoveMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

regData = getappdata(parent,'regData');
CurrentDataSet = getappdata(parent,'CurrentDataSet');

selected = get(handles.regressorGroupList,'value')-1;
selected(selected==0) =[];

regData(CurrentDataSet).regressors(selected) = [];

setappdata(parent,'CurrentRegressorGroup',selected(1)-1);
setappdata(parent,'regData',regData);

UpdateFields(hObject, eventdata, handles)

% --------------------------------------------------------------------
function splitMenu_Callback(hObject, eventdata, handles)
% hObject    handle to splitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function poolMenu_Callback(hObject, eventdata, handles)
% hObject    handle to poolMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function interactionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to interactionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

regData = getappdata(parent,'regData');
CurrentDataSet = getappdata(parent,'CurrentDataSet');

selected = get(handles.regressorGroupList,'value')-1;
selected(selected==0) =[];

regData(CurrentDataSet).regressors(end+1) = interaction(regData(CurrentDataSet).regressors(selected),'codeincr',regData(CurrentDataSet).codeincr);
regData(CurrentDataSet).codeincr = max([regData(CurrentDataSet).regressors.code]);

setappdata(parent,'CurrentRegressorGroup',length(regData(CurrentDataSet).regressors));

setappdata(parent,'regData',regData);

UpdateFields(hObject, eventdata, handles)

% --------------------------------------------------------------------
function regRemoveMenu_Callback(hObject, eventdata, handles)
% hObject    handle to regRemoveMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% parent = getappdata(handles.figure1,'parent');
% 
% regData = getappdata(parent,'regData');
% 
% selected = get(handles.regressorGroupList,'value')-1;
% selected(selected==0) =[];
% 
% regData.regressors(selected) = [];
% 
% setappdata(parent,'CurrentRegressorGroup',selected(1)-1);
% 
% Update(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function regGroupMenu_Callback(hObject, eventdata, handles)
% hObject    handle to regGroupMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function regMenu_Callback(hObject, eventdata, handles)
% hObject    handle to regMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function UpdateFields(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');

currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');
% currentRegressor = getappdata(parent,'CurrentRegressor');
CurrentDataSet = getappdata(parent,'CurrentDataSet');
regData = getappdata(parent,'regData');

if CurrentDataSet == 0
    fprintf('\nNo data set selected...')
    return
end
nreg = length(regData(CurrentDataSet).regressors);
lblnum = cellfun(@num2str, num2cell(1:nreg),'uniformoutput',false);
lbls = strcat(lblnum,{sprintf('. ')},{regData(CurrentDataSet).regressors.label});
if currentRegressorGroup > 0 &&  currentRegressorGroup <= length(regData(CurrentDataSet).regressors)
    nreg = regData(CurrentDataSet).regressors(currentRegressorGroup).Npar;
else
    nreg = 0;
end

liststrGr = cat(2,{'none selected'},lbls);
liststrR = cat(2,{'none selected'},cellfun(@num2str,mat2cell(1:nreg,1,ones(1,nreg)),'UniformOutput',0));

if get(handles.regressorGroupList,'value')> length(liststrGr)
    set(handles.regressorGroupList,'value',1)
    setappdata(parent,'CurrentRegressorGroup',0);    
end
if get(handles.regList,'value')> length(liststrR)
    set(handles.regList,'value',1)
    setappdata(parent,'CurrentRegressor',0)
end
set(handles.regressorGroupList,'string',liststrGr);
set(handles.regList,'string',liststrR);

if get(handles.displayCheck,'value')    
    DrawRegressor(hObject, eventdata, handles)
end


% --------------------------------------------------------------------
function DrawRegressor(hObject, eventdata, handles)

%Displays the value of regressors for each bin in the main axis

parent = getappdata(handles.figure1,'parent');

regData = getappdata(parent,'regData');
currentRegressorGroup =  getappdata(parent,'CurrentRegressorGroup');
currentRegressor =  getappdata(parent,'CurrentRegressor');
currentFixation =  getappdata(parent,'CurrentFixation');
currentTrial =  getappdata(parent,'CurrentTrial');
CurrentDataSet =  getappdata(parent,'CurrentDataSet');

trialData =  getappdata(parent,'trialData');
binData =  getappdata(parent,'binData');

if isempty(currentRegressor) || isequal(currentRegressor,0)  || currentRegressor > length(regData(CurrentDataSet).regressors) 
    return
end
if currentTrial > 0 && ~isempty(currentFixation) &&...
                currentFixation == 0 && ~isempty(trialData(CurrentDataSet).trials(currentTrial).fixations) 
    currentFixation = trialData(CurrentDataSet).trials(currentTrial).fixations(1);
elseif isempty(currentFixation)  || currentFixation == 0  || currentFixation > max(regData(CurrentDataSet).fixationIndex) 
    return
end
getpts = currentFixation == regData(CurrentDataSet).fixationIndex;
R = regData(CurrentDataSet).regressors(currentRegressorGroup).value(getpts,:);
if isequal(regData(CurrentDataSet).regressors(currentRegressorGroup).info.form,'sparse')
    R = unsparsify(R,'transpose');
end
tr = unique(regData(CurrentDataSet).trialIndex(getpts));

if length(tr) > 1
    error('A fixation appears to be associated with multiple trials.')
end

Zdata = R(:,currentRegressor);

if isempty(tr)
    return
end
trbincode = trialData(CurrentDataSet).trials(tr).binGroup;
bincodes = [binData.groups.code];

bindex = find(ismember(bincodes,trbincode));

Zdata = mat2cell(Zdata,cat(1,binData.groups(bindex).nbin),1);

bmfuns = getappdata(parent,'binManagerFunctions');


% crind = 0;
% bmh = getappdata(parent,'BinManager');
% bmhandles = guidata(bmh);
% set(bmhandles.binList,'value',bindex+1);

% for i = 1:length(bindex)
%     
% bgr = binData.groups(bindex(i));

%  setappdata(parent,'CurrentBinGroup',0);
 setappdata(parent,'CurrentBinGroup',bindex);

bmfuns.clear( parent );
% ph = bmfuns.draw(parent, Zdata, bindex,0);
ph = bmfuns.draw(parent, Zdata, bindex,0);

setappdata(parent,'patchHandles',ph)

% crind = crind + bgr.nbin;
% end



% % --------------------------------------------------------------------
% function  DM = designMatrix(hObject, eventdata, handles)

% --------------------------------------------------------------------
function InitializeRegressors(hObject, eventdata, handles)

%Set up the default regressors, fixation bin position, and polar
%coordinates wrt to the current and previous fixation location

parent = getappdata(handles.figure1,'parent');

crdat = getappdata(parent,'CurrentDataSet');

binData = getappdata(parent,'binData');
if isempty(binData) || isempty([binData.groups.binnums])   
    error('No bins have been defined');
end

% roiData = getappdata(parent,'roiData');
% useRois = ~(isempty(roiData) || isempty([roiData.groups.binnums]));   

trialData = getappdata(parent,'trialData');
fixData = getappdata(parent,'fixationData');
regData = getappdata(parent,'regData');



for CurrentDataSet = length(regData)+1:length(trialData)
    
    
    setappdata(parent,'CurrentDataSet',CurrentDataSet)
    if isempty(trialData) || isempty([trialData(CurrentDataSet).trials.number])   
        warning('No trials have been defined for data set %i, skipping...',CurrentDataSet);
        continue
    end
    % if  isempty([trialData(CurrentDataSet).trials.fixations])   
    %     error('No fixations are found');
    % end
    if isempty([trialData(CurrentDataSet).trials.binGroup])   
        warning('Trials in set %i have not been associated with any sampling bins, skipping...',CurrentDataSet);
        continue
    end

    if sum([trialData(CurrentDataSet).trials.nfix])==0   
        warning('Trials in set %i do not contain any fixations, skipping...',CurrentDataSet);
        continue
    end

    % X,Y positions of bin centers  for each fixation

    % posField = 'meanPos';

    % fixtrialind = zeros(nfix,1);
    % bincenters = zeros(nfix*nbin,2);

    % allfixpos = cat(1,fixData.fix.(posField));

    % fixpos = [];
    nfixs = [trialData(CurrentDataSet).trials.nfix];
    nbins = [trialData(CurrentDataSet).trials.nbin];
    nrow = sum(nfixs.*nbins);
  
%     binpos = [];
    binpos = zeros(nrow,2);
    cr = 0;
    crtinvl = 0;

    binGroupCodes = [binData.groups.code];
    nbinsfx = [];
    fixbinsovlp =[];
    logbinvolume = zeros(nrow,1);
    trialIndex = zeros(nrow,1);
    fixationIndex = zeros(nrow,1);
    
    fixnum = zeros(nrow,1);
    fixdt = zeros(nrow,1);
    fixtrt = zeros(nrow,1);
    fixbinsovlp = zeros(nrow,1);
    
    conditional_model = getappdata(parent,'ConditionalModel');
    time_intvl = getappdata(parent,'TSampIntvl');
    crtintvl = 0;
%     for i = 1:length(trialData(CurrentDataSet).trials)    
    for i = 1:length(trialData(CurrentDataSet).trials)    
        
        
        
        binIndex = ismember(binGroupCodes,trialData(CurrentDataSet).trials(i).binGroup);
        nbin = trialData(CurrentDataSet).trials(i).nbin;
        nfix = trialData(CurrentDataSet).trials(i).nfix;
        if ~isempty(fixData) && ~isequal(nfix,0)
            nbinsfx(end+(1:nfix),:) = nbin;
            fixinds = kron( trialData(CurrentDataSet).trials(i).fixations', ones(trialData(CurrentDataSet).trials(i).nbin,1) );
            
            assigninds = cr + ( 1 : nfix*nbin);
            fixnum( assigninds,: ) = fixinds-fixinds(1)+1;
            if CurrentDataSet <= length(fixData)
                fixdt( assigninds,: ) = cat(1,fixData(CurrentDataSet).fix(fixinds).dt);
                fixtrt( assigninds,: ) = cat(1,fixData(CurrentDataSet).fix(fixinds).startT) -trialData(CurrentDataSet).trials(i).startTime ;
                
%                 tbin( assigninds,: ) = ceil(fixtrt./time_intvl) + crtintvl;
                
           end            
            
            crtintvl =  crtintvl+ceil(trialData(CurrentDataSet).trials(i).stopTime./time_intvl);
            
        %     fixpos(assigninds ,:) = allfixpos(fixinds,:); % position of each fixation       
            bincent = cat(1,binData.groups(binIndex).centers);     
            fm = trialData(CurrentDataSet).trials(i).fixmat';
            fixbinsovlp( assigninds,: ) = fm(:); %cat(1,fixbinsovlp,fm(:));
            
            %Get log bin volume and normalize
            lbv = log(trialData(CurrentDataSet).trials(i).binareas);
            lbv =  lbv - repmat(mean(lbv),size(lbv,1),1);
            fxlbv = kron(ones(nfix,1),lbv(:));
%             logbinvolume = cat(1,logbinvolume,fxlbv);
            logbinvolume( assigninds,: ) =fxlbv;
            
            binpos( assigninds,: ) = repmat( bincent, nfix,1); %position of each bin center
            trialIndex( assigninds,: ) = i;
            if ~isempty(assigninds)
                fixationIndex( assigninds,: ) = fixinds; 
            end
            cr = cr + nbin*nfix;
            
        end
    end

    if ~conditional_model
       tmfun = getappdata(parent,'trialManagerFunctions');
        intvl_data = tmfun.get_time_intvl_data();
        fxind = cat(1,intvl_data.fxi)~=0;
        ntimebins = cat(1,intvl_data.fxnbins);
        expander = zeros(ntimebins'*nbinsfx,1);
        
        crind = 0;
        expi  = (1:length(fixationIndex))'; 
        for i = 1:length(fixationIndex)
            expander((1:ntimebins(i)*nbinsfx(i)) + crind,:) = expi(fixationIndex==fixationIndex(i));
            crind = crind + ntimebins(i)*nbinsfx(i);
        end
        
    else
        expander = (1:length(fixationIndex))';
    end
    
    if length(regData) < CurrentDataSet || isempty([regData(CurrentDataSet).regressors.code])
        regData(CurrentDataSet).trialIndex = trialIndex(expander);
        regData(CurrentDataSet).fixationIndex = fixationIndex(expander);
        cdi = 0;
        regData(CurrentDataSet).regressors(1) = makeregressor(full(fixbinsovlp(expander) ),'label','BinMembership','noptions',nbinsfx,'codeincr',cdi); cdi = cdi+1;
        regData(CurrentDataSet).codeincr = cdi; 
        regData(CurrentDataSet).regressors(2) = makeregressor(binpos(expander,:),'label','BinXY','noptions',nbinsfx,'codeincr',cdi);cdi = cdi+1;
         regData(CurrentDataSet).regressors(3) = makeregressor(logbinvolume,'label','BinVolume','noptions',nbinsfx,'codeincr',cdi);cdi = cdi+1;
%         regData(CurrentDataSet).regressors(3) = makeregressor(logbinvolume(expander),'label','BinVolume','noptions',nbinsfx,'codeincr',cdi);cdi = cdi+1;

    %     regData(CurrentDataSet).regressors(2) = makeregressor(getDistToBinRadial(hObject, eventdata, handles,0),'FixPolarD[0]','noptions',nbinsfx,'codeincr',regData(CurrentDataSet).codeincr);
    %     regData(CurrentDataSet).regressors(2) = makeregressor(getDistToBinRadial(hObject, eventdata, handles,-1),'label','FixPolarD[-1]','noptions',nbinsfx,'codeincr',regData(CurrentDataSet).codeincr);
    %     regData(CurrentDataSet).regressors(3) = makeregressor(getDistToBinRadial(hObject, eventdata, handles,-2),'label','FixPolarD[-2]','noptions',nbinsfx,'codeincr',regData(CurrentDataSet).codeincr);
    %     
    %     
    %     regData(CurrentDataSet).regressors(2) = makeRegressor(fixpos,...
    %                         'FixXY[0]','noptions',nbinsfx,'codeincr',regData(CurrentDataSet).codeincr);
        setappdata(parent,'regData',regData)
        regData(CurrentDataSet).regressors(4) = makeregressor(getShiftedXYFromFixation(hObject, eventdata, handles,-1,expander),...
                            'label','FixXY[-1]','noptions',nbinsfx(:)','codeincr',regData(CurrentDataSet).codeincr+2);
        regData(CurrentDataSet).regressors(5) = makeregressor(getShiftedXYFromFixation(hObject, eventdata, handles,-2,expander),...
                            'label','FixXY[-2]','noptions',nbinsfx(:)','codeincr',regData(CurrentDataSet).codeincr+3);
        regData(CurrentDataSet).regressors(6) = makeregressor(getShiftedXYFromFixation(hObject, eventdata, handles,-3,expander),...
                            'label','FixXY[-3]','noptions',nbinsfx(:)','codeincr',regData(CurrentDataSet).codeincr+4);

        regData(CurrentDataSet).regressors(7) = makeregressor(fixnum(expander,:),...
                            'label','Fixn Number','noptions',nbinsfx(:)','codeincr',regData(CurrentDataSet).codeincr+5);
        
        if conditional_model
            regData(CurrentDataSet).regressors(8) = makeregressor(fixdt,...
                            'label','Inter-Fixn T','noptions',nbinsfx(:)','codeincr',regData(CurrentDataSet).codeincr+6);
            regData(CurrentDataSet).regressors(9) = makeregressor(fixtrt,...
                            'label','Trial T','noptions',nbinsfx,'codeincr',regData(CurrentDataSet).codeincr+7);
        else
            regData(CurrentDataSet).regressors(8) = makeregressor(cat(1,intvl_data.fxtimes),...
                            'label','Inter-Fixn T','noptions',nbinsfx(:)','codeincr',regData(CurrentDataSet).codeincr+6);
            regData(CurrentDataSet).regressors(9) = makeregressor(cat(1,intvl_data.trtimes),...
                            'label','Trial T','noptions',nbinsfx,'codeincr',regData(CurrentDataSet).codeincr+7);
        end    
                        %     regData(CurrentDataSet).regressors(1) = makeRegressor(fixpos,'FixXY[0]','noptions',nbinsfx,'codeincr',regData(CurrentDataSet).codeincr);

        regData(CurrentDataSet).regressors(10) = makeregressor(cat(1,trialIndex),...
            'label','Trial Number','noptions',nbinsfx,'codeincr',regData(CurrentDataSet).codeincr+8);

        regData(CurrentDataSet).codeincr = max([regData(CurrentDataSet).regressors.code]);
%         setappdata(parent,'regData',regData)
    setappdata(parent,'regData',regData)

    end

end

setappdata(parent,'CurrentDataSet',crdat)

setappdata(parent,'regData',regData)


%----------------------------------------------------------------
function RAng = getDistToBinRadial(hObject, eventdata, handles,varargin)

% Returns the distance of each bin and angle to the (k+shift)'th fixation, where k is
% the k'th fixation event in the model. Angle is clockwise with 0
% corresponding to straight up (12 O'Clock). Units of bin position are
% assumed to have [0 0] as upper left corner

if nargin < 4
    shift = 0;
else
    shift = varargin{1};
end

expander = varargin{2};
% Returns a regressor for fixation position, shifted in sequence time by shift from the current fixation.



%Set up the default regressors

parent = getappdata(handles.figure1,'parent');

% binData = getappdata(parent,'binData');
regData = getappdata(parent,'regData');

CurrentDataSet = getappdata(parent,'CurrentDataSet');


% roiData = getappdata(parent,'roiData');
    
trialData = getappdata(parent,'trialData');
if isempty(trialData) || isempty([trialData(CurrentDataSet).trials.number])   
    error('No trials have been defined');
end
if  isempty([trialData(CurrentDataSet).trials.fixations])   
    error('No trials have been defined');
end
if isempty([trialData(CurrentDataSet).trials.binGroup])   
    error('Trials have not been associated with any sampling bins');
end



% X,Y positions of bin centers  for each fixation

posField = 'meanPos';

% fixtrialind = zeros(nfix,1);
% bincenters = zeros(nfix*nbin,2);
fixData = getappdata(parent,'fixationData');

allfixpos= cat(1,fixData(CurrentDataSet).fix.(posField));

nfixs = size(allfixpos,1);

allfixposshift = nan(size(allfixpos));

allfixposshift( max([1 , -shift+1]) : min([nfixs,nfixs-shift]),: ) =...
            allfixpos ( max([1 ,(shift+1)]) : min([nfixs,nfixs+shift]),: );
    



fixpos = [];

cr = 0;
for i = 1:length(trialData(CurrentDataSet).trials)    
    
    nbin = trialData(CurrentDataSet).trials(i).nbin;
    nfix = trialData(CurrentDataSet).trials(i).nfix;
    fixinds = kron( trialData(CurrentDataSet).trials(i).fixations', ones(trialData(CurrentDataSet).trials(i).nbin,1) );
    assigninds = cr + ( 1 : nfix*nbin);
    fixpos(assigninds ,:) = allfixposshift(fixinds,:); % position of each fixation       
   
    cr = cr + nbin*nfix;
end


% We will correct for screen aspect ratio
scrndata = getappdata(parent,'screenData');

fxnormMat = diag([1 1]./scrndata(CurrentDataSet).res(1)); %Fixation position is normalized to width

binnormMat = diag(scrndata(CurrentDataSet).res./scrndata(CurrentDataSet).res(1));

Dvec = (unsparsify( regData(CurrentDataSet).regressors(1).value,'transpose')*binnormMat - fixpos(expander,:)*fxnormMat) ;

RAng(:,1) = sqrt(sum( (Dvec).^2 ,2)); %Normalized radial distance 

RAng(:,2) = atan2(Dvec(:,1),-Dvec(:,2)); %Normalized angle



%----------------------------------------------------------------
function DFix =getShiftedXYFromFixation(hObject, eventdata, handles,varargin)

% Returns the distance of each bin and angle to the (k+shift)'th fixation, where k is
% the k'th fixation event in the model. Angle is clockwise with 0
% corresponding to straight up (12 O'Clock). Units of bin position are
% assumed to have [0 0] as upper left corner

if nargin < 4
    shift = 0;
else
    shift = varargin{1};
end

expander = varargin{2};

% Returns a regressor for fixation position, shifted in sequence time by shift from the current fixation.



%Set up the default regressors

parent = getappdata(handles.figure1,'parent');

% binData = getappdata(parent,'binData');
regData = getappdata(parent,'regData');

CurrentDataSet = getappdata(parent,'CurrentDataSet');


% roiData = getappdata(parent,'roiData');
    
trialData = getappdata(parent,'trialData');
if isempty(trialData) || isempty([trialData(CurrentDataSet).trials.number])   
    warning('No trials have been defined');
    return
end
if  isempty([trialData(CurrentDataSet).trials.fixations])   
    warning('No fixations in trial');
    return
end
if isempty([trialData(CurrentDataSet).trials.binGroup])   
    warning('Trials have not been associated with any sampling bins');
    return
end



% X,Y positions of bin centers  for each fixation

posField = 'meanPos';

% fixtrialind = zeros(nfix,1);
% bincenters = zeros(nfix*nbin,2);
fixData = getappdata(parent,'fixationData');

allfixpos= cat(1,fixData(CurrentDataSet).fix.(posField));

nfixs = size(allfixpos,1);

allfixposshift = nan(size(allfixpos));

allfixposshift( max([1 , -shift+1]) : min([nfixs,nfixs-shift]),: ) =...
            allfixpos ( max([1 ,(shift+1)]) : min([nfixs,nfixs+shift]),: );
    



fixpos = [];

cr = 0;
for i = 1:length(trialData(CurrentDataSet).trials)    
    
    nbin = trialData(CurrentDataSet).trials(i).nbin;
    nfix = trialData(CurrentDataSet).trials(i).nfix;
    fixinds = kron( trialData(CurrentDataSet).trials(i).fixations', ones(trialData(CurrentDataSet).trials(i).nbin,1) );
    assigninds = cr + ( 1 : nfix*nbin);
    fixpos(assigninds ,:) = allfixposshift(fixinds,:); % position of each fixation       
   
    cr = cr + nbin*nfix;
end


% We will correct for screen aspect ratio
scrndata = getappdata(parent,'screenData');
if CurrentDataSet> length(scrndata)
    fxnormMat = diag(scrndata(1).res.^-1); %Fixation position is normalized to width
else
    fxnormMat = diag(scrndata(CurrentDataSet).res.^-1); %Fixation position is normalized to width
end
%
% fxnormMat = diag([1 1]./scrndata.res(1)); %Fixation position is normalized to width
% 
% binnormMat = diag(scrndata.res./scrndata.res(1));

if isequal(regData(CurrentDataSet).regressors(2).info.form,'sparse')
    DFix = unsparsify( regData(CurrentDataSet).regressors(2).value,'transpose') - fixpos(expander,:)*fxnormMat ;
else
    DFix = regData(CurrentDataSet).regressors(2).value- fixpos(expander,:)*fxnormMat ;
end
% DFix(:,1) = sqrt(sum( (Dvec).^2 ,2)); %Normalized radial distance
% 
% DFix(:,2) = atan2(Dvec(:,1),-Dvec(:,2)); %Normalized angle
% 



% --- Executes on button press in displayCheck.
function displayCheck_Callback(hObject, eventdata, handles)
% hObject    handle to displayCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');
bmh = getappdata(parent,'BinManager');
if ~isempty(bmh) && ishandle(bmh) && get(handles.displayCheck,'value')
    bmhandles = guidata(bmh);
    set(bmhandles.displayCheck,'value',0)
end
% Hint: get(hObject,'Value') returns toggle state of displayCheck
 UpdateFields(hObject, eventdata, handles)



% --------------------------------------------------------------------
function modify_menu_Callback(hObject, eventdata, handles)
% hObject    handle to modify_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function PoolMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PoolMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

selected = get(handles.regressorGroupList,'value')-1;

selected(selected==0) = [];
% currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');
CurrentDataSet = getappdata(parent,'CurrentDataSet');

if isempty(selected)
    return
end
regData = getappdata(parent,'regData');

regData(CurrentDataSet).regressors(end+1) = pool(regData(CurrentDataSet).regressors(selected));
regData(CurrentDataSet).regressors(end).code = regData(CurrentDataSet).codeincr+1;
regData(CurrentDataSet).codeincr = regData(CurrentDataSet).codeincr+1;

setappdata(parent,'regData',regData)


% --------------------------------------------------------------------
function SplitMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SplitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PolynomialsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PolynomialsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SinusoidsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SinusoidsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function IneractionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to IneractionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%----------------------------------------------------
function resetRegressors(hObject,eventdata,handles)

parent = getappdata(hObject,'parent');
evh = getappdata(parent,'eyetrackerHeaderData');

for i = 1:length(evh)
    setappdata(parent,'CurrentDataSet',i)
    InitializeRegressors(hObject, eventdata, handles)
    
end
