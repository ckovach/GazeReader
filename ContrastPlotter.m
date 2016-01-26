function varargout = ContrastPlotter(varargin)
% CONTRASTPLOTTER M-file for ContrastPlotter.fig
%      CONTRASTPLOTTER, by itself, creates a new CONTRASTPLOTTER or raises the existing
%      singleton*.
%
%      H = CONTRASTPLOTTER returns the handle to a new CONTRASTPLOTTER or the handle to
%      the existing singleton*.
%
%      CONTRASTPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTRASTPLOTTER.M with the given input arguments.
%
%      CONTRASTPLOTTER('Property','Value',...) creates a new CONTRASTPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ContrastPlotter_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ContrastPlotter_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help ContrastPlotter

% Last Modified by GUIDE v2.5 17-May-2008 22:31:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ContrastPlotter_OpeningFcn, ...
                   'gui_OutputFcn',  @ContrastPlotter_OutputFcn, ...
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


% --- Executes just before ContrastPlotter is made visible.
function ContrastPlotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ContrastPlotter (see VARARGIN)



parent = varargin{1};

setappdata(handles.figure1,'parent',parent);
setappdata(handles.figure1,'ContrastPlotter',handles.figure1);

children = getappdata(parent,'children');
children(end+1) = handles.figure1;
setappdata(parent,'children',children)

plotfig = figure('visible','off');
setappdata(handles.figure1,'plotfig',plotfig);

setappdata(handles.figure1,'XaxisTerm',0);
setappdata(handles.figure1,'YaxisTerm',0);

setappdata(handles.figure1,'CurrentInteractionTerm',0);

% Choose default command line output for ContrastPlotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ContrastPlotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ContrastPlotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in plotLRR.
function plotLRR_Callback(hObject, eventdata, handles)
% hObject    handle to plotLRR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotLRR

if get(handles.plotInScene,'value')
    set([handles.LLRerror, handles.LRRoverLRRerr,handles.RRexp],'value',0);
elseif get(handles.plot1d,'value')   
    set([handles.LRRoverLRRerr,handles.RRexp],'value',0);
end


Update(hObject, eventdata, handles)

% --- Executes on button press in LLRerror.
function LLRerror_Callback(hObject, eventdata, handles)
% hObject    handle to LLRerror (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LLRerror
if get(handles.plotInScene,'value')
    set([handles.plotLRR, handles.LRRoverLRRerr,handles.RRexp],'value',0);
elseif get(handles.plot1d,'value')   
    set([handles.LRRoverLRRerr,handles.RRexp],'value',0);
end

Update(hObject, eventdata, handles)


% --- Executes on button press in RRexp.
function RRexp_Callback(hObject, eventdata, handles)
% hObject    handle to RRexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RRexp

if get(handles.plotInScene,'value')||get(handles.plot1d,'value')   
    set([handles.LLRerror,handles.plotLRR, handles.LRRoverLRRerr],'value',0);
end

Update(hObject, eventdata, handles)

% --- Executes on button press in LRRoverLRRerr.
function LRRoverLRRerr_Callback(hObject, eventdata, handles)
% hObject    handle to LRRoverLRRerr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LRRoverLRRerr

if get(handles.plotInScene,'value')||get(handles.plot1d,'value')   
    set([handles.LLRerror,handles.plotLRR,handles.RRexp],'value',0);
end

Update(hObject, eventdata, handles)

function Xval_Callback(hObject, eventdata, handles)
% hObject    handle to Xval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xval as text
%        str2double(get(hObject,'String')) returns contents of Xval as a double

Update(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function Xval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Yval_Callback(hObject, eventdata, handles)
% hObject    handle to Yval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Yval as text
%        str2double(get(hObject,'String')) returns contents of Yval as a double
Update(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function Yval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Yval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in intxnList.
function intxnList_Callback(hObject, eventdata, handles)
% hObject    handle to intxnList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns intxnList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from intxnList

selected = get(handles.intxnList,'value')-1;
% intxnindex = [0,getappdata(handles.figure1,'intxnindex')];

% selected(selected==0) = [];
% if ~isempty(intxnindex)
if ~isempty(selected)
%     setappdata(handles.figure1,'CurrentInteractionTerm',intxnindex(selected));
    setappdata(handles.figure1,'CurrentInteractionTerm',selected);
end

if isequal(get(handles.figure1,'selectiontype'),'open')
    Update(hObject, eventdata, handles)
else
    UpdateFields(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function intxnList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intxnList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function intxnValue_Callback(hObject, eventdata, handles)
% hObject    handle to intxnValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of intxnValue as text
%        str2double(get(hObject,'String')) returns contents of intxnValue as a double

parent = getappdata(handles.figure1,'parent');
CurrentDataSet = getappdata(parent,'CurrentDataSet');

selRegressorGroup =  getappdata(parent,'SelectedRegressorGroups');
regData = getappdata(parent,'regData');
RValues = getappdata(handles.figure1,'RValues');
RFunctions = getappdata(handles.figure1,'RFunctions');
% intxnRs = getappdata(handles.figure1,'intxnRs');
currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');

infos = [regData(CurrentDataSet).regressors.info];
inputTerms = {infos.functionInputCodes};
inputRegs = cellfun(@(a,b) b*ones(1,size(a,2),1),inputTerms,num2cell(1:length(infos)),'uniformoutput',false);
inputRegs = [inputRegs{:}];
inputTerms =[inputTerms{:}];
[unqterms,q,termindex] = unique(inputTerms','rows','stable');
% indx = code2ind(parent,inputTerms(1,:));

currentTerms = unique(termindex(ismember(inputRegs,selRegressorGroup )),'stable');
currentTermIndex = currentTerms(currentInteractionTerm);

RValues{currentTermIndex} = str2num(get(handles.intxnValue,'string'));

inputstr = get(handles.intxnValue,'string');
if isempty(deblank(inputstr))
    inputstr = '0';
    set(handles.intxnValue,'string',inputstr);

end
RFunctions{currentTermIndex} = str2func(['@(X) [',regexprep(inputstr,'Y','X'),']']);
setappdata(handles.figure1,'RValues',RValues);
setappdata(handles.figure1,'RFunctions',RFunctions);

Update(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function intxnValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intxnValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot1d.
function plot1d_Callback(hObject, eventdata, handles)
% hObject    handle to plot1d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot1d

if get(handles.plot1d,'value')
   set(handles.plot2d,'value',0)
   set(handles.plotInScene,'value',0)
   set(handles.Yval,'enable','off')
   set(handles.Ycheck,'enable','off')
end

Update(hObject, eventdata, handles)
% --- Executes on button press in plot2d.
function plot2d_Callback(hObject, eventdata, handles)
% hObject    handle to plot2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot2d

if get(handles.plot2d,'value')
   set(handles.plot1d,'value',0)
   set(handles.plotInScene,'value',0)
   set(handles.Yval,'enable','on')
   set(handles.Ycheck,'enable','on')
end

Update(hObject, eventdata, handles)
% --- Executes on button press in plotInScene.
function plotInScene_Callback(hObject, eventdata, handles)
% hObject    handle to plotInScene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotInScene

if get(handles.plotInScene,'value')
   set(handles.plot2d,'value',0)
   set(handles.plot1d,'value',0)
   set(handles.Yval,'enable','on')
end
Update(hObject, eventdata, handles)



%-------------------------------------------------
function Update(hObject, eventdata, handles,varargin)



parent = getappdata(handles.figure1,'parent'); %Handle for main GUI

CurrentDataSet = getappdata(parent,'CurrentDataSet');

modelData = getappdata(parent,'modelData');
regData = getappdata(parent,'regData');

currentModel = getappdata(parent,'CurrentModel');

currmodel = modelData(CurrentDataSet).models(currentModel);



%Whether the plot is 1d, 2d or on the scene plane
plottype = 1*get(handles.plot1d,'value') + 2*get(handles.plot2d,'value') + 4*get(handles.plotInScene,'value');


% currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');
selectedRegressorGroups =  getappdata(parent,'SelectedRegressorGroups');


infos = [regData(CurrentDataSet).regressors.info];
inputTerms = {infos.functionInputCodes};
inputRegs = cellfun(@(a,b) b*ones(1,size(a,2),1),inputTerms,num2cell(1:length(infos)),'uniformoutput',false);
inputRegs = [inputRegs{:}];
inputTerms =[inputTerms{:}];
[terms,q,termindex] = unique(inputTerms','rows','stable'); 

currentTerms = unique(termindex(ismember(inputRegs,selectedRegressorGroups)),'stable');


% currentTermIndex = currentTerms(currentInteractionTerm);
% inputLabels = {regData.regressors(indx).label};

setappdata(handles.figure1,'currentTerms',currentTerms);
% currentTermIndex = currentTerms(currentInteractionTerm);

modelfit = currmodel.fit(1);



RValues = getappdata(handles.figure1,'RValues');
RFunctions = getappdata(handles.figure1,'RFunctions');


rcode = [regData(CurrentDataSet).regressors.code];
if ~isequal(selectedRegressorGroups,0) && any(~ismember(rcode(selectedRegressorGroups),modelData(CurrentDataSet).models(currentModel).regressors))
    selectedRegressorGroups = 0;
    setappdata(parent,'SelectedRegressorGroups',selectedRegressorGroups);
    return
elseif selectedRegressorGroups ~= 0 
    crRegs = regData(CurrentDataSet).regressors(selectedRegressorGroups);
end

% 
% intxnRs =[];
% for i = 1:length(crRegs) 
%     a = crRegs(i).factmat(:,1);
% %     [a,unqind] = unique(a);
% %     if length(a) > 1    
% %         a(unqind) = a;
% %     end
%     intxnRs = cat(1,intxnRs,a);
% end
% 
% 
% 
% rcodes = [regData(CurrentDataSet).regressors.code];
% [a,intxnRindx] = ismember(intxnRs,rcodes);


UpdateFields(hObject, eventdata, handles)


if isempty(modelfit)
%     fprintf('\nThis model has not been fit to the data yet.')
    return
elseif isempty(currentModel) || currentModel == 0    
    return
end


plotfig = getappdata(handles.figure1,'plotfig');

if isempty(ishandle(plotfig)) || ~ishandle(plotfig)
    plotfig = figure;
    setappdata(handles.figure1,'plotfig',plotfig);
end

X = str2num(get(handles.Xval,'string'))';
Y = str2num(get(handles.Yval,'string'))';
% XY = cat(2,X,str2num(get(handles.Yval,'string'))');

activePlots = get([handles.plotLRR, handles.LLRerror, handles.LRRoverLRRerr,handles.RRexp],'value');
activePlots = [activePlots{:}]';
plotWhat = [1 2 4 8]*activePlots;

XaxisTerm = getappdata(handles.figure1,'XaxisTerm');
YaxisTerm = getappdata(handles.figure1,'YaxisTerm');


% [qq,selind] = ismember(selectedRegressorGroups,rcodes);

E = eye(modelfit.npar);

rcodes = [regData(CurrentDataSet).regressors.code];
modelreg = rcodes(ismember(rcodes,currmodel.regressors));
a = ismember(modelreg,rcodes(selectedRegressorGroups));
contrast = E(:,sum(modelfit.blockC(:,a),2)==1); 
parest = modelfit.parest; 
% infos = [regData(CurrentDataSet).regressors.info];
% inputcodes = [infos.functionInputCodes];



errmat = contrast'*modelfit.I^-1*contrast;

set(plotfig,'visible','on')
pl= [];
if length(crRegs) > 1
    crPool = pool(crRegs);
else
    crPool = crRegs;
end

ca = nan;
RV = getappdata(handles.figure1,'liststr');
RV = strcat({' '},RV(2:end));
RV(1:end-1) = strcat(RV(1:end-1),',');
RV(1) = strcat({'   '},RV(1));
pltastr = get(handles.plotArgs,'string');

if isempty(deblank(pltastr))
    plotargs = {};
else    
    plotargs = eval(pltastr);
    if ~iscell(plotargs)
        plotargs = {plotargs};
    end
end

plotvalues = getappdata(parent,'plotvalues');

plotvalues(CurrentDataSet).Rx = [];
plotvalues(CurrentDataSet).function = crPool.function;
plotvalues(CurrentDataSet).contrast = [];
plotvalues(CurrentDataSet).X = [];
plotvalues(CurrentDataSet).Y = [];
plotvalues(CurrentDataSet).Z = [];
plotvalues(CurrentDataSet).err = [];

plotvalues(CurrentDataSet).zfun = [];
plotvalues(CurrentDataSet).errfun = [];


plotvalues(CurrentDataSet).xbl = str2num(get(handles.xBaseLine,'string'));
if isempty(plotvalues(CurrentDataSet).xbl)
    blC.x = 0;
    plotvalues(CurrentDataSet).xbl=0;
else
    blC.x = 1;
end
err = [];
plotvalues(CurrentDataSet).ybl = str2num(get(handles.yBaseLine,'string'));
if isempty(plotvalues(CurrentDataSet).ybl)
    blC.y = 0;
    plotvalues(CurrentDataSet).ybl=0;
else
    blC.y = 1;
end

switch plottype
    
    case 1
        
       ca = subplot(1,1,1,'parent',plotfig);
       
       hold(ca,'on')
        
        Xi = {};
        Xbl = {}; %Baseline
        for  j = 1:length(currentTerms)
            if ismember(currentTerms(j), XaxisTerm)
                if currentTerms(j) > length(RFunctions) || ~isa(RFunctions{currentTerms(j)},'function_handle')
                    RFunctions{currentTerms(j)} = @(X) X;
                end
                Xi{j} = RFunctions{currentTerms(j)}(X);
%                 Xbl{j} = plotvalues(CurrentDataSet).xbl*blC.x;
                Xbl{j} = plotvalues(CurrentDataSet).xbl*blC.x;
            else
                Xi{j} = repmat(RValues{currentTerms(j)},size(X,1),1);
                Xbl{j} = RValues{currentTerms(j)}*blC.x;
            end
        end
        plotvalues(CurrentDataSet).Rx = Xi;
        
        regval = crPool.function(Xi{:});
        if blC.x == 1
            regval = regval-repmat(crPool.function(Xbl{:}),size(regval,1),1);
        end
        
        if bitand(plotWhat,13)
              Z = regval*(contrast'*parest);
        end
        
        if bitand(plotWhat,6)        
            err = sqrt(diag(regval*errmat*regval'));
            plotvalues(CurrentDataSet).err = err;
        end
        plotvalues(CurrentDataSet).zfun = @(parest) regval*(contrast'*parest);
%         plotvalues(CurrentDataSet).errfun = @(errmat) sqrt(diag(regval*errmat*regval'));
        plotvalues(CurrentDataSet).errfun = @(errmat)reshape(sqrt(sum( (regval*contrast'*errmat*contrast).*regval,2)),szX);
        plotvalues(CurrentDataSet).X = X;
        plotvalues(CurrentDataSet).Y = Y;
        plotvalues(CurrentDataSet).Z = Z;
        plotvalues(CurrentDataSet).contrast= contrast;
        
        if bitand(plotWhat,1) %LRR or LOR


            pl = plot(X,Z,varargin{:},plotargs{:},'parent',ca);
             title(ca(1),strcat('LRR :',RV{:}),'interpreter','none');

        end
        
        switch plotWhat
            case 2 %Error

                pl(end+1) = plot(X,err,plotargs{:},varargin{:},'parent',ca);
               title(ca(1),strcat('Err :',RV{:}),'interpreter','none');

            case 3 % LRR +/- error
        
                hold(ca,'on');
                pla = plotargs;
                pl(end+1) = plot(X,Z+err,pla{:},varargin{:},'parent',ca);
                pl(end+1) = plot(X,Z-err,pla{:},'parent',ca);
                set(pl(end-1:end),'linestyle','--');
               title(ca(1),strcat('LRR +/- err :',RV{:}),'interpreter','none');

            case 4     % Wald statistic (LRR/err)
%                 hold(ca,'off');
                pl(end+1) = plot(X,Z./err,plotargs{:},varargin{:},'parent',ca);    
                title(ca,'LRR/err')
               title(ca(1),strcat('LRR/err :',RV{:}),'interpreter','none');
            case 8       %RElative Risk
%                 hold(ca,'off');
%                 pl(end+1) = plot(X,exp(Z)./sum(exp(Z(:))),plotargs{:},varargin{:},'parent',ca);    
                pl(end+1) = plot(X,exp(Z),plotargs{:},varargin{:},'parent',ca);    
               title(ca(1),strcat('Relative Risk:',RV{:}),'interpreter','none');
            case 11 % LRR +/- error
        
                pl(end+1) = plot(X,exp(Z),plotargs{:},varargin{:},'parent',ca);    
               title(ca(1),strcat('Relative Risk +/- Err:',RV{:}),'interpreter','none');
                hold(ca,'on');
                pla = plotargs;
                pl(end+1) = plot(X,exp(Z+err),pla{:},varargin{:},'parent',ca);
                pl(end+1) = plot(X,exp(Z-err),pla{:},'parent',ca);
                set(pl(end-1:end),'linestyle','--');
               title(ca(1),strcat('LRR +/- err :',RV{:}),'interpreter','none');
        end
        
    case {2 4} % 2d plots
        
        if get(handles.polarCheck,'value') == 1
            IsRadial = true;
        else
            IsRadial = false;
        end
        
        if ~IsRadial
            [Xm,Ym] = meshgrid(X,Y);
        else
            [THm,Rm] = meshgrid(Y,X);
            Y = .25 -Y;
            
            Xm = Rm;
            Ym = THm;
%             Xm = cos(THm*2*pi).*Rm;
%             Ym = sin(THm*2*pi).*Rm;            
        end
        szXm = size(Xm);
        
        Xi = {};
%         Xbl = {};
        for  j = 1:length(currentTerms)
            if ismember(currentTerms(j),XaxisTerm) &&  ismember(currentTerms(j), YaxisTerm)              
                if currentTerms(j) > length(RFunctions) || ~isa(RFunctions{currentTerms(j)},'function_handle')
                    RFunctions{currentTerms(j)} = @(X) X;
                end
                Xi(currentTerms(j) == currentTerms) = cellfun(RFunctions{currentTerms(j)},{ Xm(:), Ym(:)},'uniformoutput',false);
%                 Xi{j} = RFunctions{currentTerms(j)}( Xm(:), Ym(:));
%                 Xbl{j} = [plotvalues(CurrentDataSet).xbl*blC.x, plotvalues(CurrentDataSet).ybl*blC.y]; 
            elseif ismember(currentTerms(j),XaxisTerm)
                if currentTerms(j) > length(RFunctions) || ~isa(RFunctions{currentTerms(j)},'function_handle')
                    RFunctions{currentTerms(j)} = @(X) X;
                end
                Xi{j} = RFunctions{currentTerms(j)}( Xm(:));
%                 Xbl{j} = plotvalues(CurrentDataSet).xbl*blC.x; 
            elseif ismember(currentTerms(j), YaxisTerm)
                if currentTerms(j) > length(RFunctions) || ~isa(RFunctions{currentTerms(j)},'function_handle')
                    RFunctions{currentTerms(j)} = @(X) X;
                end
                Xi{j} =  RFunctions{currentTerms(j)}(Ym(:));
%                 Xbl{j} = plotvalues(CurrentDataSet).ybl*blC.y;
            else
%                 Xi{j} = repmat(RValues{currentTerms(j)},numel(Xm),1);
%                 Xbl{j} = RValues{currentTerms(j)}*(blC.x | blC.y);
                Xi{j} = repmat(RValues{currentTerms(j)},numel(Xm),1);
%                 Xbl{j} = RValues{currentTerms(j)}*(blC.x | blC.y);
            end
        end
        
        plotvalues(CurrentDataSet).Rx = Xi;
        
        
        if plottype == 2    % plot on new figure
            
            nplots = sum(activePlots);
            if nplots > 1
                for k = 1:4
                	ca(k) = subplot(2,2,k,'parent',plotfig);
                end
                RV = {''};
            else
                
                ca = repmat(subplot(1,1,1,'parent',plotfig),1, 1 + 3*(nplots==1));
            end        
        elseif plottype == 4 % Plot on existing axis
%             nplots = 1;
            phandles = guidata(parent);
        	ca = repmat(phandles.axes1,1,4);
        end
        
        axis(ca,'xy');
        
        regval = crPool.function(Xi{:});
        Xbls = Xi;
        if blC.x
            xregs = ismember(currentTerms,XaxisTerm);
            Xbls(xregs) = { repmat(plotvalues(CurrentDataSet).xbl,size(regval,1),1)}; %baseline values
        end
        if blC.y
            yregs = find(ismember(currentTerms,YaxisTerm));
            if xregs == yregs
                for yr = 1:length(yregs)
                    Xbls{yregs}(:,2) = repmat(plotvalues(CurrentDataSet).ybl,size(regval,1),1); %baseline values
                end
            else
                Xbls(yregs) = { repmat(plotvalues(CurrentDataSet).ybl,size(regval,1),1)}; %baseline values
            end
        end        
        if blC.x == 1 || blC.y == 1
%              regval = regval-repmat(crPool.function(Xbl{:}),size(regval,1),1);
             regval = regval-crPool.function(Xbls{:});
         end
        
        if bitand(plotWhat,13) 
            Z = reshape(regval*(contrast'*parest),szXm);
        end    
        
        if bitand(plotWhat,6) 
            err = reshape(sqrt(sum((regval*errmat).*regval,2)),szXm);
            plotvalues(CurrentDataSet).err = err;
            
        end
        
        plotvalues(CurrentDataSet).zfun = @(parest) reshape(regval*(contrast'*parest),szXm);
        
%         plotvalues(CurrentDataSet).errfun = @(errmat) reshape(sqrt(diag(regval*contrast'*errmat*contrast*regval')),szXm);
       plotvalues(CurrentDataSet).errfun = @(errmat) reshape(sqrt(sum( (regval*contrast'*errmat*contrast).*regval,2)),szXm);


        plotvalues(CurrentDataSet).X = X;
        plotvalues(CurrentDataSet).Y = Y;
        if exist('Z','var')
            plotvalues(CurrentDataSet).Z = Z;
        end
        plotvalues(CurrentDataSet).contrast= contrast;
        
        if bitand(plotWhat , 1)  %LRR
            if IsRadial
                pcolor_radial(2*pi*Y,X,Z,ca(1));
            else         
                 pl(end+1) =imagesc(X,Y,Z,'parent',ca(1));
            end
            
             title(ca(1),strcat('LRR :',RV{:}),'interpreter','none');
               
        end
        
        if bitand(plotWhat , 2)  %1./error
            if IsRadial
                pcolor_radial(2*pi*Y,X,1./err,ca(2));
            else
                
                 pl(end+1) =imagesc(X,Y,1./err,'parent',ca(2));
            end
%             title(ca(2),'err^{-1}')
             title(ca(2),strcat('1/err :',RV{:}),'interpreter','none');

        end
        
        if bitand(plotWhat , 4)  %Wald (LRR/error)
            if IsRadial
                pcolor_radial(2*pi*Y,X,Z./err,ca(3));
            else
            
             pl(end+1) =imagesc(X,Y,Z./err,'parent',ca(3));
            end
            %             title(ca(3),'LRR/err')
             title(ca(3),strcat('LRR/err :',RV{:}),'interpreter','none');

        end
        
        if bitand(plotWhat , 8)  %Relative Risk
            if IsRadial
                pcolor_radial(2*pi*Y,X,exp(Z),ca(4));
            else
            
                 pl(end+1) =imagesc(X,Y,exp(Z),'parent',ca(4));
            end
             title(ca(4),strcat('Weight :',RV{:}),'interpreter','none');

        end
        
    otherwise
        set(plotfig,'visible','off')
end        

setappdata(parent,'plotvalues',plotvalues);

if ishandle(ca)
    axis(ca,'xy')
end
setappdata(plotfig,'plothandles',pl);


% --------------------------------------------------------------------
% function pooled_function = makePooledFunctions(grh,selReg )
% 
% 
% %Create a function which distributes inputs to the appropriate arguments
% 
% regData = getappdata(grh,'regData');
% currdat = getappdata(grh,'CurrentDataSet');
% 
% 
% 
% % Get regressor codes
% rcodes = [regData(currdat).regressors.code];
% rfunctions = {regData(currdat).regressors(selReg).function};
% factmats = {regData(currdat).regressors(selReg).factmat};
% 
% 
% terms= zeros(2,0);
% for i = 1:length(selReg)
%     [unq,q,polyterm] = unique(factmats{i}(:,1));
%      rinputs{i} = [];
%      for j = 1:length(unq)
%          for k = 1:sum(polyterm==j)
%              
%              input = find(terms(1,:) == unq(j) & terms(2,:) == k);
%              if isempty(input)
%                  terms(:,end+1) = [unq(j),k];
%                  rinputs{i}(end+1) = size(terms,2);
%              else
%                  rinputs{i}(end+1) = input;
%              end
%                  
%          end
%      end
% end             
% 
% nargs = max([rinputs{:}]);
% args = cellfun(@(n) sprintf('X%i',n),num2cell(1:nargs),'uniformoutput',false);
% trim = @(a) a(1:end-1);
% genarg = @(X)trim(sprintf('%s,',X{:}));
% funs = cellfun(@(inp,i) sprintf('rfunctions{%i}(%s)',i,genarg(args(inp))),rinputs,num2cell(1:length(rinputs)),'uniformoutput',false);
% 
% pooled_function = str2func( sprintf('@(%s) cat(2, %s )',genarg(args),genarg(funs)));
% 



% --------------------------------------------------------------------
function UpdateFields(hObject, eventdata, handles)



parent = getappdata(handles.figure1,'parent');


modelData = getappdata(parent,'modelData');
regData = getappdata(parent,'regData');

plottype = 1*get(handles.plot1d,'value') + 2*get(handles.plot2d,'value') + 4*get(handles.plotInScene,'value');


currentModel = getappdata(parent,'CurrentModel');
CurrentDataSet = getappdata(parent,'CurrentDataSet');
% currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');
selectedRegressorGroups =  getappdata(parent,'SelectedRegressorGroups');


infos = [regData(CurrentDataSet).regressors.info];
inputTerms = {infos.functionInputCodes};
inputRegs = cellfun(@(a,b) b*ones(1,size(a,2),1),inputTerms,num2cell(1:length(infos)),'uniformoutput',false);
inputRegs = [inputRegs{:}];
inputTerms =[inputTerms{:}];
[terms,q,termindex] = unique(inputTerms','rows','stable'); 

currentTerms = unique(termindex(ismember(inputRegs,selectedRegressorGroups)),'stable');



termRegs = code2ind(parent,terms(:,1));


% currentTerms = find(ismember(code2ind(parent,terms(:,1)),selectedRegressorGroups ));
% currentTermIndex = currentTerms(currentInteractionTerm);

% Create labels to display
inputLabels = {regData(CurrentDataSet).regressors(termRegs(currentTerms)).label};
for i = 1:length(currentTerms)
    crind = cumsum(termRegs(currentTerms(i)) == termRegs);    
    nterm = crind(end);
    if nterm > 1
        inputLabels{i} = sprintf('%s (X%i)',inputLabels{i},crind(currentTerms(i)));
    end
end
        


XaxisTerm = getappdata(handles.figure1,'XaxisTerm');
YaxisTerm = getappdata(handles.figure1,'YaxisTerm');

currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');
RValues = getappdata(handles.figure1,'RValues');
RFunctions = getappdata(handles.figure1,'RFunctions');

% currmodel = modelData(CurrentDataSet).models(currentModel);
% modelfit = currmodel.fit;

% 
% rcodes = [regData(CurrentDataSet).regressors.code];
% if ~isequal(selectedRegressorGroups,0) && any(~ismember(rcodes(selectedRegressorGroups),modelData(CurrentDataSet).models(currentModel).regressors))
%     selectedRegressorGroups = 0;
%     setappdata(parent,'SelectedRegressorGroups',selectedRegressorGroups);
%     return
% elseif selectedRegressorGroups~=0
%     crRegs = regData(CurrentDataSet).regressors(selectedRegressorGroups);
% end
% 
% 
% currentTerms =[];
% for i = 1:length(crRegs) 
%     a = crRegs(i).factmat(:,1);
%     [a,unqind] = unique(a);
% 
%      if length(a) > 1    
%         a(unqind) = a;
%     end
%     currentTerms = cat(1,currentTerms,a(:));
% end
% 
% 
% setappdata(handles.figure1,'currentTerms',currentTerms)
% 
% 
% 
% rcodes = [regData(CurrentDataSet).regressors.code];
% [a,currentTerms] = ismember(currentTerms,rcodes);
% 
% lbls = {regData(CurrentDataSet).regressors(currentTerms).label};

lbls = inputLabels;
sttlen = cellfun(@length,lbls);
spaces1 = cellfun(@(n) repmat(' ',1,n),mat2cell(20-2*sttlen,1,ones(1,length(sttlen))),'UniformOutput',false);

set(handles.intxnList,'Value',currentInteractionTerm+1);

RV = {};
% for i = 1:length(currentTerms)
for i = 1:length(currentTerms)
%     if ismember(currentTerms(i),XaxisTerm) && (ismember(currentTerms(i),YaxisTerm) || isequal(YaxisTerm, 0)) && plottype > 1
%         RV{i} = ' = XY';
%     elseif ismember(currentTerms(i), XaxisTerm) 
    if ismember(currentTerms(i), XaxisTerm) 
        if currentTerms(i) > length(RFunctions) || ~isa(RFunctions{currentTerms(i)},'function_handle')
                        str ='@(X)X';
                        RFunctions{currentTerms(i)} = @(X)X;
         else
                        str = char( RFunctions{currentTerms(i ) });
         end
            RV{i} =  [' = ',str(5:end)];
%         RV{i} = ' = X';

    elseif ismember(currentTerms(i), YaxisTerm) && plottype > 1
        if currentTerms(i) > length(RFunctions) || ~isa(RFunctions{currentTerms(i)},'function_handle')
                str ='@(Y)Y';
                RFunctions{currentTerms(i)} = @(Y)Y;
            
        else
                    str = char( RFunctions{currentTerms(i ) });
        end
            RV{i} =  [' = ',str(5:end)];

%         RV{i} = ' = Y';
    elseif isempty(RValues) || currentTerms(i) > length(RValues) ||  isempty(RValues{currentTerms(i)}) 
        RValues{currentTerms(i)} = zeros(1,regData(CurrentDataSet).regressors(inputRegs(currentTerms(i))).Npar);
        RV{i} = [' =', sprintf(' %0.3g',RValues{currentTerms(i)})];
    else
        RV{i} = [' =', sprintf(' %0.3g',RValues{currentTerms(i)})];
    end
end



setappdata(handles.figure1,'RValues',RValues);
setappdata(handles.figure1,'RFunctions',RFunctions);

% [liststr,intxnindex] = unique(strcat(lbls,spaces1,RV),'first');
liststr = strcat(lbls,spaces1,RV);
% setappdata(handles.figure1,'intxnindex',intxnindex);

liststr = cat(2,{'none selected'},liststr);
setappdata(handles.figure1,'liststr',liststr);

set(handles.intxnList,'string',liststr);

if currentInteractionTerm~=0  && currentInteractionTerm <= length(currentTerms)  && currentTerms(currentInteractionTerm)<=length(RValues)
%     crind = find(ismember([regData.regressors.code],currentRegressorGroup));
%     if ~isempty(crind)
    set(handles.intxnValue,'string',num2str(RValues{currentTerms(currentInteractionTerm)}));
%     end
end

if isempty(currentInteractionTerm ) %|| ~ismember(currentInteractionTerm ,intxnindex)
    currentInteractionTerm = 0;
    setappdata(handles.figure1,'CurrentInteractionTerm',currentInteractionTerm+1);
end
% set(handles.intxnList,'value',find([0,intxnindex] == currentInteractionTerm));
if currentInteractionTerm+1 <= length(liststr)
    set(handles.intxnList,'value', currentInteractionTerm+1);
else
    set(handles.intxnList,'value', length(liststr));
    currentInteractionTerm = length(liststr)-1;
end    

if currentInteractionTerm == 0
    return
end


if ismember(currentTerms(currentInteractionTerm), XaxisTerm)
    set(handles.Xcheck,'value',1);
    set(handles.Ycheck,'value',0);
%     set(handles.intxnValue,'string','X');
    if currentTerms(currentInteractionTerm) > length(RFunctions) || ~isa(RFunctions{currentTerms(currentInteractionTerm)},'function_handle')
                str ='@(X)X';
            
    else
                str = char( RFunctions{currentTerms(currentInteractionTerm ) });
    end
            set( handles.intxnValue, 'string' , str(5:end) );

elseif ismember(currentTerms(currentInteractionTerm),YaxisTerm)
    set(handles.Ycheck,'value',1);
    set(handles.Xcheck,'value',0);
%     set(handles.intxnValue,'string','Y');
    if currentTerms(currentInteractionTerm) > length(RFunctions) || ~isa(RFunctions{currentTerms(currentInteractionTerm)},'function_handle')
                str ='@(Y)Y';
            
    else
                str = char( RFunctions{currentTerms(currentInteractionTerm ) });
    end
            set( handles.intxnValue, 'string' , str(5:end) );

else   
    set(handles.Ycheck,'value',0);
    set(handles.Xcheck,'value',0);
end    
    

% --------------------------------------------------------------------
function regListContext_Callback(hObject, eventdata, handles)
% hObject    handle to regListContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SetAsPrimary_Callback(hObject, eventdata, handles)
% hObject    handle to SetAsPrimary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');


PR = get(handles.figure1,'XaxisTerm');
SR = get(handles.figure1,'YaxisTerm');

PR = unique(cat(PR,currentRegressorGroup),'stable');
SR(SR == currentRegressorGroup) = [];

set(handles.figure1,'XaxisTerm',PR);
set(handles.figure1,'YaxisTerm',SR);

Update(hObject, eventdata, handles);


% --------------------------------------------------------------------
function SetAsSecondary_Callback(hObject, eventdata, handles)
% hObject    handle to SetAsSecondary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');

currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');


PR = get(handles.figure1,'XaxisTerm');
SR = get(handles.figure1,'YaxisTerm');

SR = unique(cat(SR,currentRegressorGroup),'stable');
PR(PR == currentRegressorGroup) = [];

set(handles.figure1,'XaxisTerm',PR);
set(handles.figure1,'YaxisTerm',SR);

Update(hObject, eventdata, handles);


% --- Executes on button press in Xcheck.
function Xcheck_Callback(hObject, eventdata, handles)
% hObject    handle to Xcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Xcheck

% parent = getappdata(handles.figure1,'parent');

% currentTerms = getappdata(handles.figure1,'currentTerms');


currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');

selectedInteraction = get(handles.intxnList,'value')-1;


currentTerms = getappdata(handles.figure1,'currentTerms');

currtindx = currentTerms(selectedInteraction);

if currentInteractionTerm == 0
    return
end

PR = getappdata(handles.figure1,'XaxisTerm');
% SR = getappdata(handles.figure1,'YaxisTerm');


if get(handles.Xcheck,'value')
    set(handles.intxnValue,'string','X')
%     PR = unique(cat(2,PR,currentTerms(currentInteractionTerm)));
    PR = [PR,currtindx];
    
    RFunctions = getappdata(handles.figure1,'RFunctions');
    RFunctions{currtindx} = str2func('@(X)X');
    setappdata(handles.figure1,'RFunctions',RFunctions);
    
    setappdata(handles.figure1,'XaxisTerm',PR);
%     setappdata(handles.figure1,'YaxisTerm',SR);
%     setappdata(handles.figure1,'YaxisTerm',0);

else
    PR(PR==currtindx) = [];
%     PR(ismember(PR, currentTerms(currentInteractionTerm))) = [];
    setappdata(handles.figure1,'XaxisTerm',PR);
    
end    
    
Update(hObject, eventdata, handles);

% --- Executes on button press in Ycheck.
function Ycheck_Callback(hObject, eventdata, handles)
% hObject    handle to Ycheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ycheck

% parent = getappdata(handles.figure1,'parent');

currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');

currentTerms = getappdata(handles.figure1,'currentTerms');

currtindx = currentTerms(currentInteractionTerm);

if currentInteractionTerm == 0
    return
end

% PR = getappdata(handles.figure1,'XaxisTerm');
SR = getappdata(handles.figure1,'YaxisTerm');



if get(handles.Ycheck,'value')
    set(handles.intxnValue,'string','Y')
%     SR = cat(2,SR,currentTerms(currentInteractionTerm));
    SR = [SR,currtindx]; 
    RFunctions = getappdata(handles.figure1,'RFunctions');
    RFunctions{currentTerms(currentInteractionTerm )} = str2func('@(Y)Y');
    setappdata(handles.figure1,'RFunctions',RFunctions);

%     setappdata(handles.figure1,'XaxisTerm',PR);
%     setappdata(handles.figure1,'XaxisTerm',0);
    setappdata(handles.figure1,'YaxisTerm',SR);

else
    SR(SR==currtindx) = []; 
%     SR(ismember(SR,currentTerms(currentInteractionTerm))) = [];
    setappdata(handles.figure1,'YaxisTerm',SR);
end    
    
Update(hObject, eventdata, handles);



function plotArgs_Callback(hObject, eventdata, handles)
% hObject    handle to plotArgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plotArgs as text
%        str2double(get(hObject,'String')) returns contents of plotArgs as a double
if isequal(get(handles.figure1,'selectiontype'),'open')
    Update(hObject, eventdata, handles)
end
% --- Executes during object creation, after setting all properties.
function plotArgs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotArgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xBaseLine_Callback(hObject, eventdata, handles)
% hObject    handle to xBaseLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xBaseLine as text
%        str2double(get(hObject,'String')) returns contents of xBaseLine as a double


Update(hObject, eventdata, handles);
% --- Executes during object creation, after setting all properties.
function xBaseLine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xBaseLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yBaseLine_Callback(hObject, eventdata, handles)
% hObject    handle to yBaseLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yBaseLine as text
%        str2double(get(hObject,'String')) returns contents of yBaseLine as a double

Update(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function yBaseLine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yBaseLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function polartext_Callback(hObject, eventdata, handles)
% hObject    handle to polartext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of polartext as text
%        str2double(get(hObject,'String')) returns contents of polartext as a double


% --- Executes during object creation, after setting all properties.
function polartext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to polartext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in polarCheck.
function polarCheck_Callback(hObject, eventdata, handles)
% hObject    handle to polarCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of polarCheck

if get(handles.polarCheck,'Value') == 0;
    set(handles.text1,'string','X:');
    set(handles.text2,'string','Y:');
    set(handles.Xcheck,'string','X');
    set(handles.Ycheck,'string','Y');
else
    set(handles.text1,'string','R:');
    set(handles.text2,'string','TH:');
    set(handles.Xcheck,'string','R');
    set(handles.Ycheck,'string','TH');
end







