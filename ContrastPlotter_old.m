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

plotfig = figure('visible','off');
setappdata(handles.figure1,'plotfig',plotfig);

setappdata(handles.figure1,'PrimaryRegressor',0);
setappdata(handles.figure1,'SecondaryRegressor',0);

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

selected = get(handles.intxnList,'value');
intxnindex = [0,getappdata(handles.figure1,'intxnindex')];

% selected(selected==0) = [];
if ~isempty(intxnindex)
    setappdata(handles.figure1,'CurrentInteractionTerm',intxnindex(selected));
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
currentDataSet = getappdata(parent,'CurrentDataSet');
% currentRegressorGroup =  getappdata(parent,'CurrentRegressorGroup');
regData = getappdata(parent,'regData');
RValues = getappdata(handles.figure1,'RValues');
intxnRs = getappdata(handles.figure1,'intxnRs');
currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');
RValues{ismember([regData(currentDataSet).regressors.code],intxnRs(currentInteractionTerm))} = str2num(get(handles.intxnValue,'string'));
setappdata(handles.figure1,'RValues',RValues);

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

%Whether the plot is 1d, 2d or on the scene plane
plottype = 1*get(handles.plot1d,'value') + 2*get(handles.plot2d,'value') + 4*get(handles.plotInScene,'value');


currentModel = getappdata(parent,'CurrentModel');
currentDataSet = getappdata(parent,'CurrentDataSet');
% currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');
selectedRegressorGroups =  getappdata(parent,'SelectedRegressorGroups');


% PrimaryRegressor = getappdata(handles.figure1,'PrimaryRegressor');
% SecondaryRegressor = getappdata(handles.figure1,'SecondaryRegressor');
% currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');

RValues = getappdata(handles.figure1,'RValues');


modelData = getappdata(parent,'modelData');
regData = getappdata(parent,'regData');
currmodel = modelData(currentDataSet).models(currentModel);
modelfit = currmodel.fit(1);

% if currentRegressorGroup ~=0 && ~isempty(RValues) 
%     crind = find(ismember([regData.regressors.code],currentRegressorGroup));
%     if ~isempty(crind)
%         set(handles.intxnValue,'string',num2str(RValues{crind}));
%     end
% end

rcode = [regData(currentDataSet).regressors.code];
if ~isequal(selectedRegressorGroups,0) && any(~ismember(rcode(selectedRegressorGroups),modelData(currentDataSet).models(currentModel).regressors))
    selectedRegressorGroups = 0;
    setappdata(parent,'SelectedRegressorGroups',selectedRegressorGroups);
    return
elseif selectedRegressorGroups ~= 0 
    crRegs = regData(currentDataSet).regressors(selectedRegressorGroups);
end


intxnRs =[];
for i = 1:length(crRegs) 
    a = crRegs(i).factmat(:,1);
    [a,unqind] = unique(a);
    a(unqind) = a;
    intxnRs = cat(1,intxnRs,a);
end



rcodes = [regData(currentDataSet).regressors.code];
[a,intxnRindx] = ismember(intxnRs,rcodes);

% 
% lbls = {regData(currentDataSet).regressors(intxnRindx).label};
% sttlen = cellfun(@length,lbls);
% spaces1 = cellfun(@(n) repmat(' ',1,n),mat2cell(20-2*sttlen,1,ones(1,length(sttlen))),'UniformOutput',false);
% 
% set(handles.intxnList,'Value',currentInteractionTerm+1);

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

PrimaryRegressor = getappdata(handles.figure1,'PrimaryRegressor');
SecondaryRegressor = getappdata(handles.figure1,'SecondaryRegressor');


% [qq,selind] = ismember(selectedRegressorGroups,rcodes);

E = eye(modelfit.npar);

rcodes = [regData(currentDataSet).regressors.code];
modelreg = rcodes(ismember(rcodes,currmodel.regressors));
a = ismember(modelreg,rcodes(selectedRegressorGroups));
contrast = E(:,sum(modelfit.blockC(:,a),2)==1); 
parest = modelfit.parest; 

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

plotvalues(currentDataSet).Rx = [];
plotvalues(currentDataSet).function = crPool.function;
plotvalues(currentDataSet).contrast = [];
plotvalues(currentDataSet).X = [];
plotvalues(currentDataSet).Y = [];
plotvalues(currentDataSet).Z = [];
plotvalues(currentDataSet).err = [];


plotvalues(currentDataSet).xbl = str2num(get(handles.xBaseLine,'string'));
if isempty(plotvalues(currentDataSet).xbl)
    blC.x = 0;
    plotvalues(currentDataSet).xbl=0;
else
    blC.x = 1;
end
err = [];
plotvalues(currentDataSet).ybl = str2num(get(handles.yBaseLine,'string'));
if isempty(plotvalues(currentDataSet).ybl)
    blC.y = 0;
    plotvalues(currentDataSet).ybl=0;
else
    blC.y = 1;
end

switch plottype
    
    case 1
        
       ca = subplot(1,1,1,'parent',plotfig);
       
       hold(ca,'on')
        
        Xi = {};
        Xbl = {}; %Baseline
        for  j = 1:length(intxnRindx)
            if ismember(intxnRs(j), PrimaryRegressor)
                Xi{j} = X;
                Xbl{j} = plotvalues(currentDataSet).xbl*blC.x;
            else
                Xi{j} = repmat(RValues{intxnRindx(j)},size(X,1),1);
                Xbl{j} = RValues{intxnRindx(j)}*blC.x;
            end
        end
        plotvalues(currentDataSet).Rx = Xi;
        
        regval = crPool.function(Xi{:});
        if blC.x == 1
            regval = regval-repmat(crPool.function(Xbl{:}),size(regval,1),1);
        end
        
        if bitand(plotWhat,13)
              Z = regval*(contrast'*parest);
        end
        
        if bitand(plotWhat,6)        
            err = sqrt(diag(regval*errmat*regval'));
            plotvalues(currentDataSet).err = err;
        end
        plotvalues(currentDataSet).X = X;
        plotvalues(currentDataSet).Y = Y;
        plotvalues(currentDataSet).Z = Z;
        plotvalues(currentDataSet).contrast= contrast;
        
        if bitand(plotWhat,1) 


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
%                 strplind= cellfun(@ischar,pla);
%                 strpla = pla(strplind);
%                 if isempty(strpla), strpla = {'b'}; pla{end+1} = 'b'; end
%                 strpla(ismember(strpla,['bgcyrkw']')) = strcat(pla(ismember(strpla,['bgcyrkw']')),'--'); 
%                 pla(cellfun(@ischar,pla)) = strpla;
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
        
        Xi = {};
        Xbl = {};
        for  j = 1:length(intxnRindx)
            if ismember(intxnRs(j),PrimaryRegressor)
                Xi{j} = Xm(:);
                 Xbl{j} = plotvalues(currentDataSet).xbl*blC.x;
            elseif ismember(intxnRs(j), SecondaryRegressor)
                Xi{j} = Ym(:);
               Xbl{j} = plotvalues(currentDataSet).ybl*blC.y;
            else
%                 Xi{j} = repmat(RValues{intxnRs(j)},numel(Xm),1);
%                 Xbl{j} = RValues{intxnRs(j)}*(blC.x | blC.y);
                Xi{j} = repmat(RValues{intxnRindx(j)},numel(Xm),1);
                Xbl{j} = RValues{intxnRindx(j)}*(blC.x | blC.y);
            end
        end
        
        plotvalues(currentDataSet).Rx = Xi;
        
        
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
         if blC.x == 1 || blC.y == 1        
             regval = regval-repmat(crPool.function(Xbl{:}),size(regval,1),1);
         end
        
        if bitand(plotWhat,13) 
            Z = reshape(regval*(contrast'*parest),size(Xm));
        end    
        
        if bitand(plotWhat,6) 
            err = reshape(sqrt(diag(regval*errmat*regval')),size(Xm));
            plotvalues(currentDataSet).err = err;
            
        end

        plotvalues(currentDataSet).X = X;
        plotvalues(currentDataSet).Y = Y;
        if exist('Z','var')
            plotvalues(currentDataSet).Z = Z;
        end
        plotvalues(currentDataSet).contrast= contrast;
        
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
%              pl(end+1) =imagesc(X,Y,exp(Z)./sum(exp(Z(:))),'parent',ca(4));
            if IsRadial
                pcolor_radial(2*pi*Y,X,exp(Z),ca(4));
            else
            
                 pl(end+1) =imagesc(X,Y,exp(Z),'parent',ca(4));
            end
%             title(ca(4),'Weight')
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
function UpdateFields(hObject, eventdata, handles)



parent = getappdata(handles.figure1,'parent');

plottype = 1*get(handles.plot1d,'value') + 2*get(handles.plot2d,'value') + 4*get(handles.plotInScene,'value');


currentModel = getappdata(parent,'CurrentModel');
currentDataSet = getappdata(parent,'CurrentDataSet');
% currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');
selectedRegressorGroups =  getappdata(parent,'SelectedRegressorGroups');


PrimaryRegressor = getappdata(handles.figure1,'PrimaryRegressor');
SecondaryRegressor = getappdata(handles.figure1,'SecondaryRegressor');

currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');
RValues = getappdata(handles.figure1,'RValues');

modelData = getappdata(parent,'modelData');
regData = getappdata(parent,'regData');
% currmodel = modelData(currentDataSet).models(currentModel);
% modelfit = currmodel.fit;


rcodes = [regData(currentDataSet).regressors.code];
if ~isequal(selectedRegressorGroups,0) && any(~ismember(rcodes(selectedRegressorGroups),modelData(currentDataSet).models(currentModel).regressors))
    selectedRegressorGroups = 0;
    setappdata(parent,'SelectedRegressorGroups',selectedRegressorGroups);
    return
elseif selectedRegressorGroups~=0
    crRegs = regData(currentDataSet).regressors(selectedRegressorGroups);
end


intxnRs =[];
for i = 1:length(crRegs) 
    a = crRegs(i).factmat(:,1);
    [a,unqind] = unique(a);
    a(unqind) = a;
    intxnRs = cat(1,intxnRs,a);
end


setappdata(handles.figure1,'intxnRs',intxnRs)



rcodes = [regData(currentDataSet).regressors.code];
[a,intxnRindx] = ismember(intxnRs,rcodes);

lbls = {regData(currentDataSet).regressors(intxnRindx).label};
sttlen = cellfun(@length,lbls);
spaces1 = cellfun(@(n) repmat(' ',1,n),mat2cell(20-2*sttlen,1,ones(1,length(sttlen))),'UniformOutput',false);

set(handles.intxnList,'Value',currentInteractionTerm+1);

RV = {};
for i = 1:length(intxnRindx)
    if ismember(intxnRs(i),PrimaryRegressor) && SecondaryRegressor == 0 && plottype > 1
        RV{i} = ' = XY';
    elseif ismember(intxnRs(i), PrimaryRegressor) 
        RV{i} = ' = X';
    elseif ismember(intxnRs(i), SecondaryRegressor) && plottype > 1
        RV{i} = ' = Y';
    elseif isempty(RValues) || intxnRindx(i) > length(RValues) ||  isempty(RValues{intxnRindx(i)}) 
        RValues{intxnRindx(i)} = zeros(1,regData(currentDataSet).regressors(intxnRindx(i)).Npar);
        RV{i} = [' =', sprintf(' %0.3g',RValues{intxnRindx(i)})];
    else
        RV{i} = [' =', sprintf(' %0.3g',RValues{intxnRindx(i)})];
    end
end



setappdata(handles.figure1,'RValues',RValues);

[liststr,intxnindex] = unique(strcat(lbls,spaces1,RV),'first');
setappdata(handles.figure1,'intxnindex',intxnindex);

liststr = cat(2,{'none selected'},liststr);
setappdata(handles.figure1,'liststr',liststr);

set(handles.intxnList,'string',liststr);

if currentInteractionTerm~=0  && currentInteractionTerm <= length(intxnRs)  && intxnRs(currentInteractionTerm)<=length(RValues)
%     crind = find(ismember([regData.regressors.code],currentRegressorGroup));
%     if ~isempty(crind)
    set(handles.intxnValue,'string',num2str(RValues{intxnRs(currentInteractionTerm)}));
%     end
end

if isempty(currentInteractionTerm ) || ~ismember(currentInteractionTerm ,intxnindex)
    currentInteractionTerm = 0;
    setappdata(handles.figure1,'CurrentInteractionTerm',currentInteractionTerm);
end
set(handles.intxnList,'value',find([0,intxnindex] == currentInteractionTerm));

if currentInteractionTerm == 0
    return
end

if ismember(intxnRs(currentInteractionTerm), PrimaryRegressor)
    set(handles.Xcheck,'value',1);
    set(handles.Ycheck,'value',0);
    set(handles.intxnValue,'string','X');
elseif ismember(intxnRs(currentInteractionTerm),SecondaryRegressor)
    set(handles.Ycheck,'value',1);
    set(handles.Xcheck,'value',0);
    set(handles.intxnValue,'string','Y');
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


PR = get(handles.figure1,'PrimaryRegressor');
SR = get(handles.figure1,'SecondaryRegressor');

PR = unique(cat(PR,currentRegressorGroup));
SR(SR == currentRegressorGroup) = [];

set(handles.figure1,'PrimaryRegressor',PR);
set(handles.figure1,'SecondaryRegressor',SR);

Update(hObject, eventdata, handles);


% --------------------------------------------------------------------
function SetAsSecondary_Callback(hObject, eventdata, handles)
% hObject    handle to SetAsSecondary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');

currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');


PR = get(handles.figure1,'PrimaryRegressor');
SR = get(handles.figure1,'SecondaryRegressor');

SR = unique(cat(SR,currentRegressorGroup));
PR(PR == currentRegressorGroup) = [];

set(handles.figure1,'PrimaryRegressor',PR);
set(handles.figure1,'SecondaryRegressor',SR);

Update(hObject, eventdata, handles);


% --- Executes on button press in Xcheck.
function Xcheck_Callback(hObject, eventdata, handles)
% hObject    handle to Xcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Xcheck

% parent = getappdata(handles.figure1,'parent');

intxnRs = getappdata(handles.figure1,'intxnRs');

currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');

if currentInteractionTerm == 0
    return
end

PR = getappdata(handles.figure1,'PrimaryRegressor');
SR = getappdata(handles.figure1,'SecondaryRegressor');


if get(handles.Xcheck,'value')
    set(handles.Ycheck,'value',0)
    set(handles.intxnValue,'string','X')
    PR = unique(cat(2,PR,intxnRs(currentInteractionTerm)));
    SR(ismember(SR,intxnRs(currentInteractionTerm))) = [];
    setappdata(handles.figure1,'PrimaryRegressor',PR);
    setappdata(handles.figure1,'SecondaryRegressor',SR);

else
%     setappdata(handles.figure1,'PrimaryRegressor',0);
    PR(ismember(PR, intxnRs(currentInteractionTerm))) = [];
    setappdata(handles.figure1,'PrimaryRegressor',PR);
    
end    
    
Update(hObject, eventdata, handles);

% --- Executes on button press in Ycheck.
function Ycheck_Callback(hObject, eventdata, handles)
% hObject    handle to Ycheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ycheck


% parent = getappdata(handles.figure1,'parent');
intxnRs = getappdata(handles.figure1,'intxnRs');

currentInteractionTerm = getappdata(handles.figure1,'CurrentInteractionTerm');

% currentRegressorGroup = getappdata(parent,'CurrentRegressorGroup');

if currentInteractionTerm == 0
    return
end

PR = getappdata(handles.figure1,'PrimaryRegressor');
SR = getappdata(handles.figure1,'SecondaryRegressor');



if get(handles.Ycheck,'value')
    set(handles.Xcheck,'value',0)
    set(handles.intxnValue,'string','Y')
    SR = unique(cat(2,SR,intxnRs(currentInteractionTerm)));
    PR(ismember(PR,intxnRs(currentInteractionTerm))) = [];

    setappdata(handles.figure1,'PrimaryRegressor',PR);
    setappdata(handles.figure1,'SecondaryRegressor',SR);

else
    SR(ismember(SR,intxnRs(currentInteractionTerm))) = [];
    setappdata(handles.figure1,'SecondaryRegressor',SR);
%     setappdata(handles.figure1,'SecondaryRegressor',0);
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







