function varargout = GazeReader(varargin)
% GAZEREADER M-file for GazeReader.fig
%      GAZEREADER, by itself, creates a new GAZEREADER or raises the existing
%      singleton*.
%
%      H = GAZEREADER returns the handle to a new GAZEREADER or the handle to
%      the existing singleton*.
%
%      GAZEREADER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAZEREADER.M with the given input arguments.
%
%      GAZEREADER('Property','Value',...) creates a new GAZEREADER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GazeReader_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GazeReader_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help GazeReader

% Last Modified by GUIDE v2.5 17-Aug-2010 02:48:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GazeReader_OpeningFcn, ...
                   'gui_OutputFcn',  @GazeReader_OutputFcn, ...
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


% --- Executes just before GazeReader is made visible.
function GazeReader_OpeningFcn(hObject, eventData, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GazeReader (see VARARGIN)


CurrentDataSet = getappdata(handles.figure1,'CurrentDataSet');
if isempty(CurrentDataSet) || CurrentDataSet == 0, CurrentDataSet = 1; end


eth = getappdata(hObject,'eyetrackerHeaderData');

if CurrentDataSet <= length(eth), return, end ;  % already initialized
% Default screen properties
screenData(CurrentDataSet).res =[1024 768];
screenData(CurrentDataSet).size = 0;
screenData(CurrentDataSet).distance = 0;

%initialize data structures

hfig = handles.figure1;

% setappdata(hfig ,'buttonState',0);
setappdata(hfig ,'buttonNumber',0);
setappdata(hfig ,'roiData',[]);
setappdata(hfig ,'imageData',struct('images',makeImageStruct,'codeincr',0));
% setappdata(hfig ,'trialData',struct('trials',makeTrialData,'codeincr',0));

% setappdata(hfig ,'trialData',makeTrialData);
% setappdata(hfig ,'binData',makeBinData([]));
% setappdata(hfig ,'roiData',makeBinData([]));

setappdata(hfig ,'expEventData',makeEventData);
setappdata(hfig ,'CurrentTrial',0);
setappdata(hfig ,'CurrentFixation',0);
setappdata(hfig ,'CurrentImage',0);
setappdata(hfig ,'CurrentBinGroup',0);
setappdata(hfig ,'children',[]);
setappdata(hfig ,'activeControl','main');
setappdata(hfig ,'screenData',screenData);

%By default, assume the model is conditional on event time.
setappdata(hfig ,'ConditionalModel',true);
setappdata(hfig ,'TSampIntvl',1); % sampling interval for binnin in time



set(handles.figure1,'DeleteFcn',@Destructor)



% set(handles.figure1,'WindowButtonDownFcn',  @FigureButtonDownFcn )
% set(handles.figure1,'WindowButtonUpFcn',  @FigureButtonUpFcn )
set(handles.figure1,'units',  'pixels')

set(handles.axes2,'units',get(handles.axes1,'units'),'position',get(handles.axes1,'position'),...
            'Ydir',get(handles.axes1,'Ydir'));

axlim([2 4]) = screenData(CurrentDataSet).res;
axis(handles.axes1,axlim)
axis(handles.axes1,'manual')
% axis(handles.axes2,axlim)
axis(handles.axes2,[0 1 0 1]) %axes2 is in screen normalized coordinates
axis(handles.axes1,'manual')
% axis(handles.axes2,'manual')
% Choose default command line output for GazeReader
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GazeReader wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GazeReader_OutputFcn(hObject, eventData, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function loadData_menu_Callback(hObject, eventData, handles,varargin) %#ok
% hObject    handle to loadData_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(varargin)
    [fnames,pth] = uigetfile('*.mat','multiselect','on');
end
if length(varargin)>=1
    fnames = varargin{1};
end
if length(varargin)>=2
    pth = varargin{2};
end
    

if isnumeric(fnames)
    return
end

if ~iscell(fnames)
    fnames = {fnames};
end

D = getappdata(handles.figure1);

if isfield(D,'eyetrackerHeaderData')
    sindx0 = length(D.eyetrackerHeaderData);
else
    sindx0 = 0;
end
sindx = sindx0;


for k = 1:length(fnames)
    fname = fnames{k};
    savefname = fullfile(pth,fname);

%     matvars= whos('-file',savefname);

%     varnames = {matvars.name};


    fprintf('\nLoading %s...',savefname)
    gazeReaderData = load(savefname);
    
    if ~isfield(gazeReaderData,'gazeReaderData') && ~isfield(gazeReaderData,'eyetrackerHeaderData')              
        errordlg(sprintf('%s doesn''t\ncontain GazeReader Data',fname));
        return
    end

    
    if isfield(gazeReaderData,'gazeReaderData')
        gazeReaderData = gazeReaderData.gazeReaderData;
    end
    
    fprintf('Done.')

    fields = fieldnames(gazeReaderData); 
    dataFields = fields(~cellfun('isempty',regexp(fields,'.*Data$')));
    
    ndata(k) = length(gazeReaderData.eyetrackerHeaderData);

    
    for i = 1:length(dataFields)
        if isfield(D,(dataFields{i})) &&  ~isempty(D.(dataFields{i}))
            D.(dataFields{i})(sindx+(1:length(gazeReaderData.(dataFields{i}))) ) = gazeReaderData.(dataFields{i});
        else
            D.(dataFields{i}) = gazeReaderData.(dataFields{i});            
        end
    end
    sindx = sindx+ndata(k);
 
end

fields = fieldnames(D); 
dataFields = fields(~cellfun('isempty',regexp(fields,'.*Data$')));
for i = 1:length(dataFields)
        setappdata(handles.figure1,dataFields{i},D.(dataFields{i}));
end


if isfield(D,'trialData') %update references in trialData

    if isfield(D,'binData') && length(D.binData)>1
        [binData,codemap ] = concatBG(D.binData); 
        sindx = sindx0;
        for k = length(fnames)
            for kk = 1:ndata(k)
                newbincodes = cellfun(@(X)codemap{sindx + kk}(X),{D.trialData(sindx + kk).trials.binGroup},'uniformoutput',0);
                [D.trialData(sindx + kk).trial.binGroup] = newbincodes{:};
            end
            sindx = sindx+ndata(k);
        end
        setappdata(handles.figure1,'binData',binData);
    end
    if isfield(D,'imageData') && length(D.imageData)>1
        [imageData,codemap ] = concatIM(D.imageData); 

        sindx = sindx0;    
        for k = length(fnames)
            for kk = 1:ndata(k)
                newimcodes = cellfun(@(X)codemap{sindx + kk}(X),{D.trialData(sindx + kk).trials.image},'uniformoutput',0);
                [D.trialData(sindx + kk).trial.image] = newimcodes{:};
                sindx = sindx+ndata(k);
            end
        end
        setappdata(handles.figure1,'imageData',imageData);
    end
    setappdata(handles.figure1,'trialData',D.trialData);
end
    


if length(fnames) == 1
    setappdata(handles.figure1,'savefname',fnames{1})
end
setappdata(handles.figure1,'CurrentDataSet',0)
viewDataSetsMenu_Callback(hObject, eventData, handles);

% --------------------------------------------------------------------
function loadImage_menu_Callback(hObject, eventData, handles) %#ok<DEFNU>
% hObject    handle to loadImage_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s
% --------------------------------------------------------------------
function File_Callback(hObject, eventData, handles)
% hObject    handle to File (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function samplingBins_menu_Callback(hObject, eventData, handles)
% hObject    handle to samplingBins_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tools_menu_Callback(hObject, eventData, handles)
% hObject    handle to tools_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventData, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function grid_menu_Callback(hObject, eventData, handles)
% hObject    handle to grid_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



if ~isappdata(handles.figure1,'gridManager') || isempty(getappdata(handles.figure1,'gridManager')) || ~ishandle(getappdata(handles.figure1,'gridManager')) 
    h = GridManager(handles.figure1);
%     [h,gridManagerFunctions] = GridManager(handles.figure1);
%     setappdata(handles.figure1,'gridManager',h);
%     setappdata(handles.figure1,'gridManagerFunctions',gridManagerFunctions);
    children = getappdata(handles.figure1,'children');
    children(end+1) = h;
    setappdata(handles.figure1,'children',children)
else
    h = getappdata(handles.figure1,'gridManager');
    figure(h)
     setappdata(handles.figure1,'activeControl','gridManager');

end

 figurePositionManager(h, handles,0)

% --------------------------------------------------------------------
function samplingBins_Rectangles_menu_Callback(hObject, eventData, handles)
% hObject    handle to samplingBins_Rectangles_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function polygon_Menu_Callback(hObject, eventData, handles)
% hObject    handle to polygon_Menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Ellipses_Callback(hObject, eventData, handles)
% hObject    handle to Ellipses (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function modelFitting_menu_Callback(hObject, eventData, handles)
% hObject    handle to modelFitting_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function RoiManagerMenu_Callback(hObject, eventData, handles)
% hObject    handle to RoiManagerMenu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ch = getappdata(handles.figure1,'children');

ch(end+1) = RoiManager(handles.figure1);

figurePositionManager(ch(end),handles);

setappdata(handles.figure1,'children',ch);



% --------------------------------------------------------------------
function Regressors_menu_Callback(hObject, eventData, handles)
% hObject    handle to Regressors_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ch = getappdata(handles.figure1,'children');
 
ch(end+1) = RegManager(handles.figure1);


setappdata(handles.figure1,'children',ch);

figure(ch(end))

% --------------------------------------------------------------------
function new_Menu_Callback(hObject, eventData, handles)
% hObject    handle to new_Menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_menu_Callback(hObject, eventData, handles)
% hObject    handle to save_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

savefname = getappdata(handles.figure1,'savefname');


if isempty(savefname)
    saveAs_menu_Callback(hObject, eventData, handles)
else
    D = getappdata(handles.figure1);     %#ok<NASGU>
    
    fields = fieldnames(D);
    dataFields = fields(~cellfun('isempty',regexp(fields,'.*Data$')));

    for i = 1:length(dataFields)
        gazeReaderData.(dataFields{i}) = D.(dataFields{i});
    end
%     save(savefname,'gazeReaderData','-v7.3');
    save(savefname,'-struct','gazeReaderData','-v7.3');
end

% --------------------------------------------------------------------
function saveAs_menu_Callback(hObject, eventData, handles,fname,pth)
% hObject    handle to saveAs_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles)
    handles = guidata(hObject);
end

if nargin < 5
    pth = cd;
end

if nargin < 4
    [fname,pth] = uiputfile('*.mat');
end

if isnumeric(fname)
   return
end


savefname = fullfile(pth,fname);    
setappdata(handles.figure1,'savefname',savefname);
save_menu_Callback(hObject, eventData, handles)
    



% --------------------------------------------------------------------
function trial_Menu_Callback(hObject, eventData, handles)
% hObject    handle to trial_Menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Trial_Manager_Callback(hObject, eventData, handles)
    trialManager_menu_Callback(hObject, eventData, handles)
    
% --------------------------------------------------------------------
function trialManager_menu_Callback(hObject, eventData, handles)
% hObject    handle to trialManager_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ch = getappdata(handles.figure1,'children');

ch(end+1) = TrialManager(handles.figure1);

setappdata(handles.figure1,'children',ch);
figure(ch(end))

% --------------------------------------------------------------------
function import_menu_Callback(hObject, eventData, handles)
% hObject    handle to import_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function importEyetrackerData_menu_Callback(hObject, eventData, handles,filenames,fpath)
% hObject    handle to importEyetrackerData_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


supported_file_types = {'edf','fix','eyd'};

dataFileDir = getappdata(handles.figure1,'dataFileDir');

if nargin < 4 || isempty(filenames)
    [filenames,fpath] = uigetfile({sprintf('*.%s;',supported_file_types{:})},'Select a Data File', dataFileDir ,'multiselect','on');
elseif nargin < 5
    fpath = cd;
end

if isnumeric(filenames)
    return
end

setappdata(handles.figure1,'dataFileDir',fpath)

eyetrackerHeaderData = getappdata(handles.figure1,'eyetrackerHeaderData');
trialData = getappdata(handles.figure1,'trialData');
expEventData= getappdata(handles.figure1,'expEventData');
currentDataSet = 1 + length(eyetrackerHeaderData);
setappdata(handles.figure1,'CurrentDataSet',currentDataSet);


etfixdata=getappdata(handles.figure1,'fixationData');
if isempty(etfixdata)
    clear etfixdata;
end
    
etrawdata=getappdata(handles.figure1,'rawGazeData');
if isempty(etrawdata)
    clear etrawdata;
end


if ~iscell(filenames)
    filenames = {filenames};
end

screenData = getappdata(handles.figure1,'screenData');
for i = 1:length(filenames)
    filename = filenames{i};

    [a,b,ext] = fileparts(filename);

    switch lower(ext)

        case '.edf'
            [fixdata,rawdata] = readEDF(fullfile(fpath,filename),'subtract_first_xdat_time');        
        case '.fix'
            fixdata = readFIX(fullfile(fpath,filename));        
            try
                rawdata = ReadEYD(fullfile(fpath,regexprep(filename,'\.fix$','\.eyd'))); 
            catch
                rawdata = struct([]);
                warning(sprintf('Unable to find a .eyd file with the same name as this .fix file.\nOnly fixation data is loaded.'))
            end
        case '.eyd'
            rawdata = ReadEYD(fullfile(fpath,filename));        
            try
                fixdata = readFIX(fullfile(fpath,regexprep(filename,'\.eyd$','\.fix'))); 
            catch
                fixdata = struct([]);
                warning(sprintf('Unable to find a .fix file for this raw .eyd file.\nOnly raw data is loaded.'))
            end
        otherwise
                error('Unsupported file type.');
    end

%     if isempty(eyetrackerData)
%          eyetrackerData = struct('raw',struct('seg',[]),'fix',struct('seg',[]));
%     end
    
    nseg = max([length(fixdata.seg),length(rawdata.seg)]);

   rf = fieldnames(rawdata);    
   rf(strcmp(rf,'seg')) = [];
   fxf = fieldnames(rawdata);    
   fxf(strcmp(rf,'seg')) = [];
   
    for j = 1:nseg
%         eyetrackerData(currentDataSet + i + j - 2) = struct('raw',[],'fix',[]);
        if j <= length(rawdata.seg)
%             eyetrackerData(currentDataSet + i + j - 2).raw = rawdata.seg(j);
          etrawdata(currentDataSet + i + j - 2) = rawdata.seg(j);       
          for f = 1:length(rf)
              eyetrackerHeaderData(currentDataSet + i + j - 2).info.rawheader.(rf{f}) = rawdata.(rf{f});
          end          
        end
        if j <= length(fixdata.seg)
%             eyetrackerData(currentDataSet + i + j - 2).fix = fixdata.seg(j);
            etfixdata(currentDataSet + i + j - 2) = fixdata.seg(j);
            for f = 1:length(fxf)
              eyetrackerHeaderData(currentDataSet + i + j - 2).info.fixheader.(fxf{f}) = rawdata.(fxf{f});
            end          
        end
        
        xdatTs = [fixdata.seg(j).xdat.startT];
        xdatcodenum = [fixdata.seg(j).xdat.id];
        type = zeros(1,length(xdatTs));
        xdatlabel = {fixdata.seg(j).xdat.code};

        eyetrackerHeaderData(currentDataSet + i + j - 2).filename = filename;
        eyetrackerHeaderData(currentDataSet + i + j - 2).label= sprintf('Data %i',currentDataSet + i-1);
        eyetrackerHeaderData(currentDataSet + i + j - 2).filetype = ext(2:end);
        eyetrackerHeaderData(currentDataSet + i + j - 2).info= [];
        eyetrackerHeaderData(currentDataSet + i + j - 2).units= fixdata.units;
        
        if currentDataSet + i + j - 2 > length(expEventData)
            expEventData(currentDataSet + i + j - 2) = makeEventData('time',xdatTs,'type',type,'xdatcode',xdatcodenum,'label',xdatlabel);
        else
            expEventData(currentDataSet + i + j - 2) = makeEventData(expEventData(currentDataSet + i + j - 2),'time',xdatTs,'type',type,'xdatcode',xdatcodenum,'label',xdatlabel);
        end
        expEventData(currentDataSet + i + j - 2).xdat = fixdata.seg(j).xdat;
        
        if currentDataSet + i + j - 2 > length(trialData)
            if isempty(trialData)
                trialData = makeTrialData([]);
            else
                trialData(currentDataSet + i + j - 2) = makeTrialData([]);
            end
        end
        
%         expEventData(currentDataSet).events = newevents.events;    
%         expEventData(currentDataSet).codeincr = length(xdatTs);

        screenData(currentDataSet + i + j - 2) = screenData(end);

    end
    
    
end

% % setappdata(handles.figure1,'eyetrackerData',eyetrackerData);
% setappdata(handles.figure1,'fixationData',fixdata);
% setappdata(handles.figure1,'rawGazeData',rawdata);
setappdata(handles.figure1,'fixationData',etfixdata);
setappdata(handles.figure1,'rawGazeData',etrawdata);
setappdata(handles.figure1,'eyetrackerHeaderData',eyetrackerHeaderData);
setappdata(handles.figure1,'expEventData',expEventData);
setappdata(handles.figure1,'trialData',trialData);
setappdata(handles.figure1,'screenData',screenData);
setappdata(handles.figure1,'CurrentTrial',0);
setappdata(handles.figure1,'CurrentFixation',0);

viewDataSetsMenu_Callback(hObject,eventData,handles);



% --------------------------------------------------------------------
function importMatFile_menu_Callback(hObject, eventData, handles, varargin)
% hObject    handle to importMatFile_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% fixvnames = {'fixdata','FIX'}; %variable names that contain the fixation and raw data structures
% rawvnames = {'rawdata','RAW'};

if isempty(handles)
    handles = guidata(hObject);
end

dataFileDir = getappdata(handles.figure1,'dataFileDir');


if nargin > 4 && ~isempty(varargin{2})
    fpath= varargin{2};
    dataFileDir = fpath;
else
        fpath = [cd,filesep];
end

if nargin > 3 && ~isempty(varargin{1})
    filenames = varargin{1};
    if ~iscell(filenames)
        filenames = {filenames};
    end
else
    [filenames,fpath] = uigetfile('mat','Select Data File(s)', dataFileDir ,'multiselect','on');
end

    

if isnumeric(filenames)
    return
end

etfixdata=getappdata(handles.figure1,'fixationData');
if isempty(etfixdata)
    clear etfixdata;
end
    
etrawdata=getappdata(handles.figure1,'rawGazeData');
if isempty(etrawdata)
    clear etrawdata;
end

setappdata(handles.figure1,'dataFileDir',fpath)

eyetrackerHeaderData = getappdata(handles.figure1,'eyetrackerHeaderData');
trialData = getappdata(handles.figure1,'trialData');
expEventData= getappdata(handles.figure1,'expEventData');
currentDataSet = 1 + length(eyetrackerHeaderData);
setappdata(handles.figure1,'CurrentDataSet',currentDataSet);
screenData = getappdata(handles.figure1,'screenData');

if ~iscell(filenames)
    filenames = {filenames};
end

for i = 1:length(filenames)
    filename = filenames{i};

    load(fullfile(fpath,filename));
    
%     if exist('fixdata','var') 
%         ;
    if exist('FIX','var') 
       fixdata = FIX;
       clear FIX
%     else        
%         return
    end
%     if exist('rawdata','var') 
%         ;
    if exist('RAW','var') 
       rawdata = RAW;
       if isfield(RAW,'segData')
           rawdata.seg = RAW.segData;
       end
       clear RAW
%     else
%         return
    end
    
    if ~exist('rawdata','var') && ~exist('fixdata','var')
        warning('The data file %s does not contain the variables I''m looking for.',filename)
        continue
    end
    
    nseg = max([length(fixdata.seg),length(rawdata.seg)]);

   rf = fieldnames(rawdata);    
   rf(strcmp(rf,'seg')) = [];
   fxf = fieldnames(rawdata);    
   fxf(strcmp(rf,'seg')) = [];
   
    for j = 1:nseg
%         eyetrackerData(currentDataSet + i + j - 2) = struct('raw',[],'fix',[]);
        if j <= length(rawdata.seg)
%             eyetrackerData(currentDataSet + i + j - 2).raw = rawdata.seg(j);
          etrawdata(currentDataSet + i + j - 2) = rawdata.seg(j);       
          for f = 1:length(rf)
              eyetrackerHeaderData(currentDataSet + i + j - 2).info.rawheader.(rf{f}) = rawdata.(rf{f});
          end          
        end
        if j <= length(fixdata.seg)
%             eyetrackerData(currentDataSet + i + j - 2).fix = fixdata.seg(j);
            etfixdata(currentDataSet + i + j - 2) = fixdata.seg(j);
            for f = 1:length(fxf)
              eyetrackerHeaderData(currentDataSet + i + j - 2).info.fixheader.(fxf{f}) = rawdata.(fxf{f});
            end          
        end
        
        xdatTs = [fixdata.seg(j).xdat.startT];
        xdatcodenum = [fixdata.seg(j).xdat.id];
        type = zeros(1,length(xdatTs));
        xdatlabel = {fixdata.seg(j).xdat.code};

        eyetrackerHeaderData(currentDataSet + i + j - 2).filename = filename;
        eyetrackerHeaderData(currentDataSet + i + j - 2).label= sprintf('Data %i',currentDataSet + i-1);
        eyetrackerHeaderData(currentDataSet + i + j - 2).filetype = 'mat';
        eyetrackerHeaderData(currentDataSet + i + j - 2).info= [];
        eyetrackerHeaderData(currentDataSet + i + j - 2).units= fixdata.units;
        screenData(currentDataSet + i + j - 2)= screenData(end);
        
        if currentDataSet + i + j - 2 > length(expEventData)
            expEventData(currentDataSet + i + j - 2) = makeEventData('time',xdatTs,'type',type,'xdatcode',xdatcodenum,'label',xdatlabel);
        else
            expEventData(currentDataSet + i + j - 2) = makeEventData(expEventData(currentDataSet + i + j - 2),'time',xdatTs,'type',type,'xdatcode',xdatcodenum,'label',xdatlabel);
        end
        expEventData(currentDataSet + i + j - 2).xdat = fixdata.seg(j).xdat;
        
        if isempty(trialData)
            trialData = makeTrialData([]);            
        elseif currentDataSet + i + j - 2 > length(trialData)
            trialData(currentDataSet + i + j - 2) = makeTrialData([]);
        end
        
%         expEventData(currentDataSet).events = newevents.events;    
%         expEventData(currentDataSet).codeincr = length(xdatTs);


    end
end

% % setappdata(handles.figure1,'eyetrackerData',eyetrackerData);
% setappdata(handles.figure1,'fixationData',fixdata);
% setappdata(handles.figure1,'rawGazeData',rawdata);
setappdata(handles.figure1,'fixationData',etfixdata);
setappdata(handles.figure1,'rawGazeData',etrawdata);
setappdata(handles.figure1,'eyetrackerHeaderData',eyetrackerHeaderData);
setappdata(handles.figure1,'expEventData',expEventData);
setappdata(handles.figure1,'trialData',trialData);
setappdata(handles.figure1,'screenData',screenData);

viewDataSetsMenu_Callback(hObject,eventData,handles);

% --------------------------------------------------------------------
function Text_menu_Callback(hObject, eventData, handles)
% hObject    handle to Text_menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% %------------------------------------------------
%   
% function FigureButtonDownFcn(hObject,eventData)


%------------------------------------------------
  
function gridManagerIsAcctive_ButtonDnFcn(hObject,eventData,handles)

% Callback for the main figure window when the Grid Manager is Active


handles = guidata(hObject);

axhandles = cat(1,handles.axes1,get(handles.axes1,'children'),...
                  handles.axes2,get(handles.axes2,'children'));              
                   
if ~ismember(gco,axhandles)
    return
end


if ~isequal(get(handles.figure1,'selectiontype'),'normal') %Only left button single click will activate resizing
    return
end
currentBinGroup= getappdata(handles.figure1,'CurrentBinGroup');
binData = getappdata(handles.figure1,'binData');
% imageData = getappdata(handles.figure1,'imageData');
% screenData = getappdata(handles.figure1,'screenData');
if currentBinGroup==0 || ~strcmp(binData.groups(currentBinGroup).type,'grid')
    activefigure(handles.figure1)
    return
end

axlim = axis(handles.axes2);

if ishandle(getappdata(handles.figure1,'gridManager'))
    gmhandles = guidata(getappdata(handles.figure1,'gridManager'));
else
    gmhandles = guidata(gridManager(handles.figure1));
    setappdata(handles.figure1,'GridManager',gmhandles)
    binData = getappdata(handles.figure1,'binData');

end

if isempty(binData.groups(currentBinGroup).pos )
%     newbin = makeBinData({[0 1 0 1],[15 15]},'type','grid','label','New Grid');
%     binData.groups(currentBinGroup) = newbin.groups;
    binData = makeBinData(binData,{[0 1 0 1],[15 15]},'type','grid','label','New Grid');
end            

if currentBinGroup ~= 0    
    
    if ~isempty(binData.groups(currentBinGroup).inputData)
        newpos = repositionBox(binData.groups(currentBinGroup).inputData, handles.axes2);    
    else
        newpos = [0 1 0 1];
    end
    nbins = str2num([get(gmhandles.Nx,'string'),'  ',get(gmhandles.Ny,'string')]); %#ok<ST2NM>
    binData.groups(currentBinGroup).inputData = {newpos./axlim([2 2 4 4])};
    binData.groups(currentBinGroup).inputData{2} = nbins;
    binData.groups(currentBinGroup).pos = grid2rect(binData.groups(currentBinGroup).inputData);
end


setappdata(handles.figure1,'binData',binData);
gridfuns = getappdata(handles.figure1,'gridManagerFunctions');

% gridfuns.updateFields(hObject,eventData,guidata(getappdata(handles.figure1,'gridManager')));
gridfuns.updateFields(handles.figure1);

if ishandle(getappdata(handles.figure1,'BinManager'))
    setappdata(handles.figure1,'CurrentBin',0)
    binfuns = getappdata(handles.figure1,'binManagerFunctions');
%     binfuns.update(hObject,eventData,guidata(getappdata(handles.figure1,'BinManager')))
    binfuns.update();
else
%     gridfuns.updateBinData(hObject,eventData,guidata(getappdata(handles.figure1,'gridManager')));
%     gridfuns.updatePlot(hObject,eventData,guidata(getappdata(handles.figure1,'gridManager')));
    gridfuns.updatePlot();
end
activefigure(handles.figure1);

%------------------------------------------------

function binManagerIsAcctive_ButtonDnFcn(hObject,eventData, handles)

% Callback for the main figure window when the BinManager is Active
bm = getappdata(handles.figure1,'BinManager');
bmh = guidata(bm);
if get(bmh.edit_bins_check,'value')==1
    MoveBins(hObject, eventData,handles);
end

activefigure(handles.figure1);

%------------------------------------------------

function roiManagerIsAcctive_ButtonDnFcn(hObject,eventData, handles)

% Callback for the main figure window when the BinManager is Active

MoveRois(hObject, eventData,handles);

activefigure(handles.figure1);

%------------------------------------------------
  
% function FigureButtonUpFcn(hObject,eventData)
% 
% setappdata(hObject,'buttonState',0)
% setappdata(hObject,'buttonNumber',0)


% --------------------------------------------------------------------
function Image_Manager_Callback(hObject, eventData, handles)

         Images_Menu_Callback(hObject, eventData, handles)
% --------------------------------------------------------------------
function Images_Menu_Callback(hObject, eventData, handles)
% hObject    handle to Images_Menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~isappdata(handles.figure1,'ImageManager') || isempty(getappdata(handles.figure1,'ImageManager')) || ~ishandle(getappdata(handles.figure1,'ImageManager')) 
    h = ImageManager(handles.figure1);
%     setappdata(handles.figure1,'ImageManager',h);
% %     setappdata(handles.figure1,'ImageManagerFunctions',ImageManagerFunctions);
%     children = getappdata(handles.figure1,'children');
%     children(end+1) = h;
%     setappdata(handles.figure1,'children',children)
else
    h = getappdata(handles.figure1,'ImageManager');
    figure(h)
    setappdata(handles.figure1,'activeControl','ImageManager');

end

 figurePositionManager(h, handles,0)

% --------------------------------------------------------------------
function   ImageManagerIsAcctive_ButtonDnFcn(hObject,eventData,handles)

% 
%

                   
axhandles = cat(1,handles.axes1,get(handles.axes1,'children'),...
                  handles.axes2,get(handles.axes2,'children'));              
                   
if ~ismember(gco,axhandles)
    return
end

imhandles = guidata(getappdata(handles.figure1,'ImageManager'));


if get(imhandles.scaleCheckBox,'value') %No repositioning if the image is scaled to screen
    return
end

imfuns =  getappdata(handles.figure1,'ImageManagerFunctions');

if ~isequal(get(handles.figure1,'selectiontype'),'normal')
    return
end
currentImage = getappdata(handles.figure1,'CurrentImage');
imageData = getappdata(handles.figure1,'imageData');

if currentImage ~= 0    
    imhandles = guidata(getappdata(handles.figure1,'ImageManager'));
    newpos = repositionBox(imageData.images(currentImage).position, handles.axes2);    
    dpd = abs(diff(imageData.images(currentImage).position  - newpos));
    if ~get(imhandles.pixelScaleCheckBox,'value') || (dpd(1) < 1e-10 && dpd(3) < 1e-10 )
        imageData.images(currentImage).position = newpos;
    end
%     imageData.images(currentImage).position = newpos./imageData.images(currentImage).screenres([1 1 2 2]);
    
end
 
setappdata(handles.figure1,'imageData',imageData);
% imfuns.UpdateImage(hObject,eventData,guidata(getappdata(handles.figure1,'ImageManager')));
imfuns.UpdateImage();
    
% --------------------------------------------------------------------
function newaxrect = repositionBox(axrectnormin, ax)

if iscell(axrectnormin)
    axrectnormin = axrectnormin{1};
end

selectionTol = 10;

% A function for repositioning rectangular areas using rbbox and dragrect-
% converting from figure corrdinates to axis coordinates

currpt = get(ax, 'CurrentPoint');
currpt = currpt(1,1:2);

set(ax,'units', get(get(ax,'parent'),'units') )

axpos = get(ax,'position'); %Position of axis in figure

% screenData = getappdata(handles.figure1,'screenData');

minmax = cat(1,min(axrectnormin,[],1),max(axrectnormin,[],1)); %use only extreme corners for groups of rectangles
axrectnorm = minmax([1 4 5 8]);

% axlim = screenData.res;
axlim = axis(ax);

axrect= axrectnorm.*axlim([2 2 4 4]);
% axrect = axrect([1 2 4 3]); %reverse x and y to be consistent with xfig

xyfig = [axpos(1) axpos(2); axpos(1)+axpos(3), axpos(2);...
          axpos(1:2)+axpos(3:4); axpos(1), axpos(2)+axpos(4)];   

xyax = [axlim([1 4]) ;axlim([2 4]);  axlim([2 3]); axlim([1 3])];

xyax(:,end+1)= 1;
 xyfig(:,end+1) = 1;

% ax2fig = (xyax'*xyax)^-1*(xyax'*xyfig(:,1:2));
% fig2ax = (xyfig'*xyfig)^-1*(xyfig'*xyax(:,1:2));
ax2fig = (xyax'*xyax)^-1*(xyax'*xyfig);
% fig2ax = (xyfig'*xyfig)^-1*(xyfig'*xyax);

figrectmat =   [axrect([1 4])  1;axrect([2 3]) 1]*ax2fig;

figrect = figrectmat([1 2 3 4]);

figpt =  [currpt 1]*ax2fig ;

corners = [figrect([1 3]); figrect([2 3]); figrect([2 4]); figrect([1 4])];

d = sqrt(sum( (corners - repmat(figpt(:,1:2),4,1)).^2,2));

[mind,minc] = min(d);

cornerSelect = 0;
if mind < selectionTol
    cornerSelect = minc;
end

figrectInCornerWidthFormat = figrect([1 3 2 4]) - [0,0,figrect([1 3])];

if cornerSelect ==0
    
    newfigrect = dragrect(figrectInCornerWidthFormat);
    
else
    newfigrect = rbbox(figrectInCornerWidthFormat,corners(mod(cornerSelect+1,4)+1,:));
end

newfigmat = [newfigrect([1 2]);newfigrect([1 2])+[newfigrect(3) 0 ] ;...
             newfigrect([1 2])+newfigrect(3:4); newfigrect([1 2]) + [0 newfigrect(4) ]];
newfigmat(:,end+1) = 1;
newaxmat = newfigmat*ax2fig^-1;
newaxrect = newaxmat([1 2 8 5]);            

% newaxrect = newaxrect([1 2 4 3]); %Back into original format


% --------------------------------------------------------------------
function MoveBins(hObject, eventData,handles)

%Updates plots of the selected bin while the button is down
% and adjusts position based on figure current point
selectionTol = 10;

binData = getappdata(handles.figure1,'binData');
% screenData = getappdata(handles.figure1,'screenData');
currentBinGroup = getappdata(handles.figure1,'CurrentBinGroup');
% currentBin = getappdata(handles.figure1,'CurrentBin');


if isempty(currentBinGroup) ||currentBinGroup == 0 
    return
elseif strcmp(binData.groups(currentBinGroup).type,'grid')
    gridManagerIsAcctive_ButtonDnFcn(hObject,eventData,handles)
    return
end

set(handles.axes1,'units', get(handles.figure1,'units'))
set(handles.axes2,'ydir',get(handles.axes1,'ydir'),'units', get(handles.figure1,'units'))

axpos = get(handles.axes2,'position'); %in figure coordinates
% axlim = axis(handles.axes2)./screenData.res([1 1 2 2]);
axlim = [0 1 0 1]; %axes2 fixed in normalized coordinates

xyfig = [axpos(1) axpos(2); axpos(1)+axpos(3), axpos(2);...
          axpos(1:2)+axpos(3:4); axpos(1), axpos(2)+axpos(4)];   
xyax = [axlim([1 4]) ;axlim([2 4]);  axlim([2 3]); axlim([1 3])];
xyax(:,end+1)= 1;
xyfig(:,end+1) = 1;
ax2fig = (xyax'*xyax)^-1*(xyax'*xyfig);

currpt = get(handles.axes2, 'CurrentPoint');
% currpt = currpt(1,1:2)./screenData.res;
currpt = currpt(1,1:2);
% figpt =  [currpt 1]*ax2fig ;

% d = sqrt(sum( (corners - repmat(figpt(:,1:2),4,1)).^2,2));
% [mind,minc] = min(d);

type = binData.groups(currentBinGroup).type;

if ~isempty(binData.groups(currentBinGroup).pos)
    
    insidebins =  binData.groups(currentBinGroup).isinside( binData.groups(currentBinGroup), currpt,selectionTol/sqrt(abs(det(ax2fig))));

elseif isempty(binData.groups(currentBinGroup).pos) 
    code = binData.groups(currentBinGroup).code;
    switch type
        case 'grid'
            return
        case 'rect'
            newbindata = makeBinData(currpt([1 1 2 2]),'type',type,'label',binData.groups(currentBinGroup).label,'code',code );
        case 'ellipse'
            newbindata = makeBinData([currpt 0],'type',type,'label',binData.groups(currentBinGroup).label,'code',code );
        case 'poly'
            newbindata = makeBinData(currpt,'type',type,'label',binData.groups(currentBinGroup).label,'code',code );            
        case {'tri','simplex'}
            newbindata = makeBinData({[1 2 3], currpt([1 1 1],:)},'type',type,'label',binData.groups(currentBinGroup).label,'code',code );            
    end
    
    binData.groups(currentBinGroup) = newbindata.groups;
    insidebins = 1;
% else
%     
%     insidebins =  binData.groups(currentBinGroup).isinside( binData.groups(currentBinGroup), currpt);
%     if ~any(insidebins)
%         return
%     end
end

if ~any(insidebins)
    pos = binData.groups(currentBinGroup).pos;
    code = binData.groups(currentBinGroup).code;

    switch type
        case 'grid'
            return
        case 'rect'
            newbindata = makeBinData(cat(1,pos,currpt([1 1 2 2])),'type',type,'label',binData.groups(currentBinGroup).label,'code',code);
        case 'ellipse'
            newbindata = makeBinData(cat(1,pos,[currpt 0 0 0 ]),'type',type,'label',binData.groups(currentBinGroup).label,'code',code);
        case {'tri','simplex'}
            newpos = pos;
            nv = size(newpos.vert,1);
            [srt,srti] = sort(sum((pos.vert-repmat(currpt,nv,1)).^2,2));
            newpos.vert = cat(1,newpos.vert,currpt);
            newpos.tri = cat(1,newpos.tri,[srti(1:2)',nv+1]);
            newbindata = makeBinData(newpos,'type',type,'label',binData.groups(currentBinGroup).label,'code',code);
        case 'poly'
            newbindata = makeBinData(cat(1,pos,currpt),'type',type,'label',binData.groups(currentBinGroup).label,'code',code);            
    end
    
    binData.groups(currentBinGroup) = newbindata.groups;
    insidebins(size(binData.groups(currentBinGroup).pos,1)) = true;

end

% type = binData.groups(currentBinGroup).type;
% 
switch type
    case 'poly' 
    pos = binData.groups(currentBinGroup).pos;
    currentBin = 0;
%     dists = sqrt(sum((pos- repmat(currpt,size(pos,1),1)).^2,2));
%     [mind, currentVertex ] = min(dists);
    case {'tri','simplex'}
        pos = binData.groups(currentBinGroup).pos;
        currentBin = 0;
        
    otherwise
        insidecenters = binData.groups(currentBinGroup).centers(insidebins,:);
        insideindex = find(insidebins);
        dists = sqrt(sum(insidecenters - repmat(currpt,size(insidecenters,1),1).^2,2));
        [mind, ind] = min(dists);
        currentBin = insideindex(ind);
        setappdata(handles.figure1,'CurrentBin',currentBin);
        pos = binData.groups(currentBinGroup).pos(currentBin,:);
end    
%Find the closest bin in the group



startfpt = get(handles.figure1,'CurrentPoint');
% startfpt = getappdata(handles.figure1,'CurrentFigPoint');

axes(handles.axes2)
hold(handles.axes2,'on')
% plotbin = currentBin;
plotbin = 0;
phs = binData.groups(currentBinGroup).plot(binData.groups(currentBinGroup),plotbin,'g--','parent',handles.axes2);
% axes(handles.axes2)
bmfuns = getappdata(handles.figure1,'binManagerFunctions');

switch type
    
    case 'grid' %already handled above
  %%%%%%%%%%%%%%%
    case 'rect'
        
        corners = reshape(pos([1 4 2 4 2 3 1 3]),2,4)';
        corners(:,end+1) =1;
        
        fcorners = corners*ax2fig;
        currfpt = fcorners;
        
        [minfd,minfc] = min( sqrt(sum((fcorners(:,1:2) - repmat(startfpt ,4,1)).^2,2)));
        
        modindex = mod(minfc - 1 + ([0 1 3] ),4)+1; 
            
        while getappdata(handles.figure1,'ButtonDown')
            pause(.05)

            newfpt = get(handles.figure1,'CurrentPoint');
%             newfpt = getappdata(handles.figure1,'CurrentFigPoint');
            D = newfpt - startfpt ;
            
            if minfd < selectionTol
                if mod(minfc,2) == 1
                    currfpt( modindex,1:2) = cat(1,D, flipud(diag(D))) + fcorners( modindex,1:2) ;
                else
                    currfpt( modindex,1:2) = cat(1,D, diag(D)) + fcorners(modindex,1:2);
                end            
            else
                currfpt(:,1:2) = fcorners(:,1:2) + repmat(D,4,1);
            end
            corners = currfpt/ax2fig;
            rect = corners([1 2 7 6]);
            rect = [min(rect(1:2)) max(rect(1:2)),...
                    min(rect(3:4)) max(rect(3:4))];
                
            binData.groups(currentBinGroup).pos(currentBin,:) = rect;
            delete(phs)
            phs = binData.groups(currentBinGroup).plot(binData.groups(currentBinGroup),plotbin,'g--','parent',handles.axes2);
            axis(handles.axes2,'off')
            
            drawnow    
        end
%         0
%     bmfuns.update();
%%%%%%%%%%%%%%%%
    case 'ellipse'
        
        %transformation from normalized ellipse coordinates to axis
        %coordinates
        trmat = diag([pos(3:4),1])*[cos(pos(5)), -sin(pos(5)),   0  
                                     sin(pos(5))  cos(pos(5))   0 
                                    pos(1)        pos(2)    1];        
%         itrmat = diag(pos(3:4).^-1)*[cos(pos(5)), sin(pos(5)); -sin(pos(5))  cos(pos(5))];
%          dcent = sqrt(sum( (( [currpt-pos(1:2),0] )/trmat).^2,2)); 
        if ~det(trmat) == 0 
             ptnrm =  [currpt,1]/trmat; %point in ellipse normalized units
            
             drimnrm = [ptnrm(1:2) - ptnrm(1:2)./sqrt(sum(ptnrm.^2)-1),0]; %radial distance vector to the edge of the bin in normalized units
             drimax = drimnrm*trmat; %radial distance to the edge of the bin in axis  units
              dcent = sqrt(sum(drimax.^2));
%              trmat
%              sqrt(sum(ptnrm(1:2).^2))
%              return
        else
             dcent = 0;
             trmat = eye(3);
             trmat(3,1:2) = pos(1:2);
%              ptnrm = [0 0 1];
 
        end
        currpos = pos;
        while getappdata(handles.figure1,'ButtonDown')
            
            newfpt = get(handles.figure1,'CurrentPoint');
%             newfpt = getappdata(handles.figure1,'CurrentFigPoint');
            Df = newfpt - startfpt ;
%             posf = [pos(1:2),1]*ax2fig;
            newpt  = [newfpt,1]/ax2fig;
            newptnrm  = newpt/trmat; % new point in ellipse normalized coordinates
            if ~any(newptnrm(1:2))
                newptnrm   = [1 0 1]/(trmat*ax2fig);
            end
          
%             if abs(dcent - pos(3)) < selectionTol/sqrt(abs(det(ax2fig)))
            if dcent < selectionTol/sqrt(abs(det(ax2fig)))
                if strcmp(get(handles.figure1,'selectiontype'),'normal')
                    sax = cat(2,diag(newptnrm(1:2)),[0 0]')*trmat;%new semi axes in axis coordinates                    
%                     slengths = sqrt(sum(sax.^2,2)-1);
                    slengths = abs(sum(sax,2));
%                     currpos = pos + [0 0 slengths' 0 ]
                    currpos(3:4) = slengths' ;
                elseif strcmp(get(handles.figure1,'selectiontype'),'extend')
                    
                    th = atan2(newptnrm(2),-newptnrm(1));
                    currpos = pos + [0 0 0 0 th ];
                else
                    return
                end
            elseif strcmp(get(handles.figure1,'selectiontype'),'normal')

                currpos(1:3) = pos(1:3) + [Df,0]/ax2fig;
                
            else 
                return
            end
            binData.groups(currentBinGroup).pos(currentBin,:) = currpos;
            delete(phs)
            phs = binData.groups(currentBinGroup).plot(binData.groups(currentBinGroup),plotbin,'g--','parent',handles.axes2);
%             axis(handles.axes2,axlim);
            axis(handles.axes2,'off');
%             gca
%             getappdata(handles.figure1,'CurrentBin')
            drawnow    
%             axis(handles.axes2)

        end                
%         bmfuns.update();
%%%%%%%%%%%%%
    case 'poly'

        
        pos(:,end+1) =1;
        fvertices = pos*ax2fig;
        stpos = pos;
        [minfd,minfc] = min(sqrt(sum((fvertices - repmat([startfpt,1] ,size(pos,1),1)).^2)));
        
        while getappdata(handles.figure1,'ButtonDown')

            newfpt = get(handles.figure1,'CurrentPoint');
%             newfpt = getappdata(handles.figure1,'CurrentFigPoint');

            D = newfpt - startfpt ;
            
            if minfd < selectionTol          
                newpt = [newfpt,1]/ax2fig;
                pos(minfc,:) = newpt;
            else
                pos = repmat([D,0]/ax2fig,size(pos,1),1)+stpos;
            end
            
            binData.groups(currentBinGroup).pos = pos(:,1:2);
            delete(phs)
            phs = binData.groups(currentBinGroup).plot(binData.groups(currentBinGroup),plotbin,'g--','parent',handles.axes2);
            axis(handles.axes2,'off')
            drawnow    
        end
%         bmfuns.update();
    
%%%%%%%%%%%%%%
    case {'tri','simplex'}
%         
%         pos(:,end+1) =1;
%         fvertices = pos*ax2fig;
%         stpos = pos;
%         [minfd,minfc] = min(sqrt(sum((fvertices - repmat([startfpt,1] ,size(pos,1),1)).^2)));
%    
        vert = pos.vert;
        vert(:,end+1) = 1;
        fvertices = vert*ax2fig;
        
        nv = size(pos.vert,1);
        [minfd,minfc] = sort(sum((pos.vert-repmat(startfpt,nv,1).^2)));
        

        while getappdata(handles.figure1,'ButtonDown')

            newfpt = get(handles.figure1,'CurrentPoint');
%             newfpt = getappdata(handles.figure1,'CurrentFigPoint');

            D = newfpt - startfpt ;
            
            if minfd(1) < selectionTol          
                newpt = [newfpt,1]/ax2fig;
                pos.vert(minfc(1),:) = newpt(:,1:2);
            else
                pos.vert = repmat([D,0]/ax2fig,size(vert,1),1)+stpos;
            end
            
            binData.groups(currentBinGroup).pos = pos;
            delete(phs)
            phs = binData.groups(currentBinGroup).plot(binData.groups(currentBinGroup),plotbin,'g--','parent',handles.axes2);
            axis(handles.axes2,'off')
            drawnow    
        end
    
end

setappdata(handles.figure1,'binData',binData)
bmfuns.update();
delete(phs(ishandle(phs)))

% --------------------------------------------------------------------
function MoveRois(hObject, eventData,handles)

%Updates plots of the selected roi while the button is down
% and adjusts position based on figure current point
selectionTol = 10;

roiData = getappdata(handles.figure1,'roiData');
% screenData = getappdata(handles.figure1,'screenData');
currentRoiGroup = getappdata(handles.figure1,'CurrentRoiGroup');
% currentRoi = getappdata(handles.figure1,'CurrentRoi');


if currentRoiGroup == 0
    return
elseif strcmp(roiData.groups(currentRoiGroup).type,'grid')
    gridManagerIsAcctive_ButtonDnFcn(hObject,eventData,handles)
    return
end

set(handles.axes1,'units', get(handles.figure1,'units'))
set(handles.axes2,'ydir',get(handles.axes1,'ydir'),'units', get(handles.figure1,'units'))

axpos = get(handles.axes2,'position'); %in figure coordinates
% axlim = axis(handles.axes2)./screenData.res([1 1 2 2]);
axlim = [0 1 0 1]; %axes2 fixed in normalized coordinates

xyfig = [axpos(1) axpos(2); axpos(1)+axpos(3), axpos(2);...
          axpos(1:2)+axpos(3:4); axpos(1), axpos(2)+axpos(4)];   
xyax = [axlim([1 4]) ;axlim([2 4]);  axlim([2 3]); axlim([1 3])];
xyax(:,end+1)= 1;
xyfig(:,end+1) = 1;
ax2fig = (xyax'*xyax)^-1*(xyax'*xyfig);

currpt = get(handles.axes2, 'CurrentPoint');
% currpt = currpt(1,1:2)./screenData.res;
currpt = currpt(1,1:2);
% figpt =  [currpt 1]*ax2fig ;

% d = sqrt(sum( (corners - repmat(figpt(:,1:2),4,1)).^2,2));
% [mind,minc] = min(d);

type = roiData.groups(currentRoiGroup).type;

if ~isempty(roiData.groups(currentRoiGroup).pos)
    
    insiderois =  roiData.groups(currentRoiGroup).isinside( roiData.groups(currentRoiGroup), currpt,selectionTol/sqrt(abs(det(ax2fig))));

elseif isempty(roiData.groups(currentRoiGroup).pos) 
    code = roiData.groups(currentRoiGroup).code;
    switch type
        case 'grid'
            return
        case 'rect'
            newroidata = makeBinData(currpt([1 1 2 2]),'type',type,'label',roiData.groups(currentRoiGroup).label,'code',code );
        case 'ellipse'
            newroidata = makeBinData([currpt 0],'type',type,'label',roiData.groups(currentRoiGroup).label,'code',code );
        case 'poly'
            newroidata = makeBinData(currpt,'type',type,'label',roiData.groups(currentRoiGroup).label,'code',code );            
    end
    
    roiData.groups(currentRoiGroup) = newroidata.groups;
    insiderois = 1;
% else
%     
%     insiderois =  roiData.groups(currentRoiGroup).isinside( roiData.groups(currentRoiGroup), currpt);
%     if ~any(insiderois)
%         return
%     end
end

if ~any(insiderois)
    pos = roiData.groups(currentRoiGroup).pos;
    code = roiData.groups(currentRoiGroup).code;

    switch type
        case 'grid'
            return
        case 'rect'
            newroidata = makeBinData(cat(1,pos,currpt([1 1 2 2])),'type',type,'label',roiData.groups(currentRoiGroup).label,'code',code);
        case 'ellipse'
            newroidata = makeBinData(cat(1,pos,[currpt 0 0 0 ]),'type',type,'label',roiData.groups(currentRoiGroup).label,'code',code);
        case 'poly'
            newroidata = makeBinData(cat(1,pos,currpt),'type',type,'label',roiData.groups(currentRoiGroup).label,'code',code);            
    end
    
    roiData.groups(currentRoiGroup) = newroidata.groups;
    insiderois(size(roiData.groups(currentRoiGroup).pos,1)) = true;

end

% type = roiData.groups(currentRoiGroup).type;
% 
if ~strcmp(type,'poly') 
    insidecenters = roiData.groups(currentRoiGroup).centers(insiderois,:);
    insideindex = find(insiderois);
    dists = sqrt(sum(insidecenters - repmat(currpt,size(insidecenters,1),1).^2,2));
    [mind, ind] = min(dists);
    currentRoi = insideindex(ind);
    setappdata(handles.figure1,'CurrentRoi',currentRoi);
    pos = roiData.groups(currentRoiGroup).pos(currentRoi,:);
else
    pos = roiData.groups(currentRoiGroup).pos;
    currentRoi = 0;
%     dists = sqrt(sum((pos- repmat(currpt,size(pos,1),1)).^2,2));
%     [mind, currentVertex ] = min(dists);
end    
%Find the closest roi in the group



startfpt = get(handles.figure1,'CurrentPoint');
% startfpt = getappdata(handles.figure1,'CurrentFigPoint');

axes(handles.axes2)
hold(handles.axes2,'on')
% plotroi = currentRoi;
plotroi = 0;
phs = roiData.groups(currentRoiGroup).plot(roiData.groups(currentRoiGroup),plotroi,'g--','parent',handles.axes2);
% axes(handles.axes2)
roifuns = getappdata(handles.figure1,'roiManagerFunctions');

switch type
    
    case 'grid' %already handled above
    case 'rect'
        
        corners = reshape(pos([1 4 2 4 2 3 1 3]),2,4)';
        corners(:,end+1) =1;
        
        fcorners = corners*ax2fig;
        currfpt = fcorners;
        
        [minfd,minfc] = min( sqrt(sum((fcorners(:,1:2) - repmat(startfpt ,4,1)).^2,2)));
        
        modindex = mod(minfc - 1 + ([0 1 3] ),4)+1; 
            
        while getappdata(handles.figure1,'ButtonDown')
            pause(.05)

            newfpt = get(handles.figure1,'CurrentPoint');
%             newfpt = getappdata(handles.figure1,'CurrentFigPoint');
            D = newfpt - startfpt ;
            
            if minfd < selectionTol
                if mod(minfc,2) == 1
                    currfpt( modindex,1:2) = cat(1,D, flipud(diag(D))) + fcorners( modindex,1:2) ;
                else
                    currfpt( modindex,1:2) = cat(1,D, diag(D)) + fcorners(modindex,1:2);
                end            
            else
                currfpt(:,1:2) = fcorners(:,1:2) + repmat(D,4,1);
            end
            corners = currfpt/ax2fig;
            rect = corners([1 2 7 6]);
            rect = [min(rect(1:2)) max(rect(1:2)),...
                    min(rect(3:4)) max(rect(3:4))];
                
            roiData.groups(currentRoiGroup).pos(currentRoi,:) = rect;
            delete(phs)
            phs = roiData.groups(currentRoiGroup).plot(roiData.groups(currentRoiGroup),plotroi,'g--','parent',handles.axes2);
            axis(handles.axes2,'off')
            
            drawnow    
        end
%         0
%     roifuns.update();

    case 'ellipse'
        
        %transformation from normalized ellipse coordinates to axis
        %coordinates
        trmat = diag([pos(3:4),1])*[cos(pos(5)), -sin(pos(5)),   0  
                                     sin(pos(5))  cos(pos(5))   0 
                                    pos(1)        pos(2)    1];        
%         itrmat = diag(pos(3:4).^-1)*[cos(pos(5)), sin(pos(5)); -sin(pos(5))  cos(pos(5))];
%          dcent = sqrt(sum( (( [currpt-pos(1:2),0] )/trmat).^2,2)); 
        if ~det(trmat) == 0 
             ptnrm =  [currpt,1]/trmat; %point in ellipse normalized units
            
             drimnrm = [ptnrm(1:2) - ptnrm(1:2)./sqrt(sum(ptnrm.^2)-1),0]; %radial distance vector to the edge of the roi in normalized units
             drimax = drimnrm*trmat; %radial distance to the edge of the roi in axis  units
              dcent = sqrt(sum(drimax.^2));
%              trmat
%              sqrt(sum(ptnrm(1:2).^2))
%              return
        else
             dcent = 0;
             trmat = eye(3);
             trmat(3,1:2) = pos(1:2);
%              ptnrm = [0 0 1];
 
        end
        currpos = pos;
        while getappdata(handles.figure1,'ButtonDown')
            
            newfpt = get(handles.figure1,'CurrentPoint');
%             newfpt = getappdata(handles.figure1,'CurrentFigPoint');
            Df = newfpt - startfpt ;
%             posf = [pos(1:2),1]*ax2fig;
            newpt  = [newfpt,1]/ax2fig;
            newptnrm  = newpt/trmat; % new point in ellipse normalized coordinates
            if ~any(newptnrm(1:2))
                newptnrm   = [1 0 1]/(trmat*ax2fig);
            end
          
%             if abs(dcent - pos(3)) < selectionTol/sqrt(abs(det(ax2fig)))
            if dcent < selectionTol/sqrt(abs(det(ax2fig)))
                if strcmp(get(handles.figure1,'selectiontype'),'normal')
                    sax = cat(2,diag(newptnrm(1:2)),[0 0]')*trmat;%new semi axes in axis coordinates                    
%                     slengths = sqrt(sum(sax.^2,2)-1);
                    slengths = abs(sum(sax,2));
%                     currpos = pos + [0 0 slengths' 0 ]
                    currpos(3:4) = slengths' ;
                elseif strcmp(get(handles.figure1,'selectiontype'),'extend')
                    
                    th = atan2(newptnrm(2),-newptnrm(1));
                    currpos = pos + [0 0 0 0 th ];
                else
                    return
                end
            elseif strcmp(get(handles.figure1,'selectiontype'),'normal')

                currpos(1:3) = pos(1:3) + [Df,0]/ax2fig;
                
            else 
                return
            end
            roiData.groups(currentRoiGroup).pos(currentRoi,:) = currpos;
            delete(phs)
            phs = roiData.groups(currentRoiGroup).plot(roiData.groups(currentRoiGroup),plotroi,'g--','parent',handles.axes2);
%             axis(handles.axes2,axlim);
            axis(handles.axes2,'off');
%             gca
%             getappdata(handles.figure1,'CurrentRoi')
            drawnow    
%             axis(handles.axes2)

        end                
%         roifuns.update();

    case 'poly'

        
        pos(:,end+1) =1;
        fvertices = pos*ax2fig;
        stpos = pos;
        [minfd,minfc] = min(sqrt(sum((fvertices - repmat([startfpt,1] ,size(pos,1),1)).^2)));
        
        while getappdata(handles.figure1,'ButtonDown')

            newfpt = get(handles.figure1,'CurrentPoint');
%             newfpt = getappdata(handles.figure1,'CurrentFigPoint');

            D = newfpt - startfpt ;
            
            if minfd < selectionTol          
                newpt = [newfpt,1]/ax2fig;
                pos(minfc,:) = newpt;
            else
                pos = repmat([D,0]/ax2fig,size(pos,1),1)+stpos;
            end
            
            roiData.groups(currentRoiGroup).pos = pos(:,1:2);
            delete(phs)
            phs = roiData.groups(currentRoiGroup).plot(roiData.groups(currentRoiGroup),plotroi,'g--','parent',handles.axes2);
            axis(handles.axes2,'off')
            drawnow    
        end
%         roifuns.update();

end

setappdata(handles.figure1,'roiData',roiData)
roifuns.update();
delete(phs(ishandle(phs)))



% --------------------------------------------------------------------
function figurePositionManager(figh, handles,vert)

%Decides where to put new windows

if nargin < 3
    vert = 0;
end
rim = 55;

set(figh,'units',get(handles.figure1,'units'));

figpos = get(figh,'position');
mainfigpos = get(handles.figure1,'position');
scrsize = get(0,'screensize');
% pos = [0 1]

newpos = figpos;
if vert
    if scrsize(4) - sum(mainfigpos([2 4])) <= mainfigpos(2)
        newpos(1:2) = mainfigpos(1:2)+[mainfigpos(3)-figpos(3),-figpos(4)-rim] ;
    else
        newpos(1:2) = mainfigpos(1:2)+[mainfigpos(3)-figpos(3),mainfigpos(4)+rim] ;
    end

else
    if scrsize(3) - sum(mainfigpos([1 3])) <= mainfigpos(1)
        newpos(1:2) = mainfigpos(1:2)+[-figpos(3)-rim,mainfigpos(4)-figpos(4)] ;
    else
        newpos(1:2) = mainfigpos(1:2)+[mainfigpos(3)+rim,mainfigpos(4)-figpos(4)] ;
    end
end
set(figh,'position',newpos)
        


% --------------------------------------------------------------------
function Destructor(hObject, eventData)
% hObject    handle to Images_Menu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

children = getappdata(hObject,'children');

for i = 1:length(children)
    if ishandle(children(i)), delete(children(i));end
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventData, handles)
% hObject    handle to axes1 (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventData, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventData, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventData, handles)
% hObject    handle to Untitled_12 (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function binManagerMenu_Callback(hObject, eventData, handles)
% hObject    handle to binManagerMenu (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

binmgr = BinManager(handles.figure1);
children = getappdata(handles.figure1,'children');
children(end+1) = binmgr;
setappdata(handles.figure1,'children',children)
setappdata(handles.figure1,'BinManager',binmgr)

figurePositionManager(binmgr, handles,0)

figure(binmgr)

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventData, handles)
% hObject    handle to figure1 (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%What to do when a button is pressed over the window
handles = guidata(hObject);
setappdata(handles.figure1,'ButtonDown',1)
switch getappdata(hObject,'activeControl')
    
    case 'gridManager'
         gridManagerIsAcctive_ButtonDnFcn(hObject,eventData, handles);
    case 'ImageManager'
          ImageManagerIsAcctive_ButtonDnFcn(hObject,eventData, handles);
    case 'binManager'
          binManagerIsAcctive_ButtonDnFcn(hObject,eventData, handles);
    case 'roiManager'        
          roiManagerIsAcctive_ButtonDnFcn(hObject,eventData, handles);        
    case 'main'
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventData, handles)
% hObject    handle to figure1 (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setappdata(handles.figure1,'ButtonDown',0);


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventData, handles)
% hObject    handle to figure1 (see GCBO)
% eventData  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% setappdata(handles.figure1,'CurrentFigPoint',get(handles.figure1,'CurrentPoint'));
get(handles.figure1,'CurrentPoint');


% --------------------------------------------------------------------
function viewDataSetsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to viewDataSetsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ch = getappdata(handles.figure1,'children');

ch(end+1) = DataSetWindow(handles.figure1);

setappdata(handles.figure1,'children',ch);
setappdata(handles.figure1,'DataSetWindow',ch(end));


% --------------------------------------------------------------------
function EventManager_Callback(hObject, eventdata, handles)
% hObject    handle to EventManager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles)
    handles = guidata(hObject);
end
    
ch = getappdata(handles.figure1,'children');
% ch(end+1) = TrialBuilder(parent);
ch(end+1) = EventManager(handles.figure1);
setappdata(handles.figure1,'children',ch);


% --------------------------------------------------------------------
function models_menu_Callback(hObject, eventdata, handles)
% hObject    handle to models_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ch = getappdata(handles.figure1,'children');

ch(end+1) = ModelManager(handles.figure1);

setappdata(handles.figure1,'children',ch);




% --------------------------------------------------------------------
function DataPlotter_Callback(hObject, eventdata, handles)
% hObject    handle to DataPlotter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ch = getappdata(handles.figure1,'children');

ch(end+1) = DataDisplayManager(handles.figure1);

setappdata(handles.figure1,'children',ch);

figure(ch(end))


% --------------------------------------------------------------------
function plotContrast_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plotContrast_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ch = getappdata(handles.figure1,'children');
ch(end+1) = ContrastPlotter(handles.figure1);
setappdata(handles.figure1,'children',ch)

figure(ch(end))
% --------------------------------------------------------------------
function plot1dcontrast_Callback(hObject, eventdata, handles)
% hObject    handle to plot1dcontrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
return

function plot2Dcontrast_menu_Callback(hObject, eventdata, handles)

return 


% --------------------------------------------------------------------
