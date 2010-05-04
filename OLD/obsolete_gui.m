function varargout = obsolete_gui(varargin)
% OBSOLETE_GUI M-file for obsolete_gui.fig
%      OBSOLETE_GUI, by itself, creates a new OBSOLETE_GUI or raises the existing
%      singleton*.
%
%      H = OBSOLETE_GUI returns the handle to a new OBSOLETE_GUI or the handle to
%      the existing singleton*.
%
%      OBSOLETE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OBSOLETE_GUI.M with the given input arguments.
%
%      OBSOLETE_GUI('Property','Value',...) creates a new OBSOLETE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before obsolete_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to obsolete_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help obsolete_gui

% Last Modified by GUIDE v2.5 10-Dec-2007 16:42:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @obsolete_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @obsolete_gui_OutputFcn, ...
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


% --- Executes just before obsolete_gui is made visible.
function obsolete_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to obsolete_gui (see VARARGIN)

setappdata(hObject,'buttonState',0)
setappdata(hObject,'buttonNumber',0)
setappdata('roiData',[])
setappdata('trialData',[])
    
set(handles.figure1,'WindowButtonDownFcn',  @FigureButtonDownFcn )
set(handles.figure1,'WindowButtonUpFcn',  @FigureButtonUpFcn )

% Choose default command line output for obsolete_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes obsolete_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = obsolete_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function xrange1_Callback(hObject, eventdata, handles)
% hObject    handle to xrange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xrange1 as text
%        str2double(get(hObject,'String')) returns contents of xrange1 as a double


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

% Hints: get(hObject,'String') returns contents of xrange2 as text
%        str2double(get(hObject,'String')) returns contents of xrange2 as a double


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

% Hints: get(hObject,'String') returns contents of yrange1 as text
%        str2double(get(hObject,'String')) returns contents of yrange1 as a double


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

% Hints: get(hObject,'String') returns contents of yrange2 as text
%        str2double(get(hObject,'String')) returns contents of yrange2 as a double


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

% Hints: get(hObject,'String') returns contents of Nx as text
%        str2double(get(hObject,'String')) returns contents of Nx as a double


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



function Ny_Callback(hObject, eventdata, handles)
% hObject    handle to Ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ny as text
%        str2double(get(hObject,'String')) returns contents of Ny as a double


% --- Executes during object creation, after setting all properties.
function Ny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function loadImage_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to loadImage_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imfmts = imformats;
imexts = strcat('*.',[imfmts.ext]);

[img,pth] = uigetfile( {sprintf('%s;',imexts{:}),'Image Files'})


trialData(ROIdata.currentIndex).imagefile = img;
trialData(ROIdata.currentIndex).imagepath = pth;

axes(handles.axes1);

setappdata(hObject,'trialData',trialData)

guiData(hObject);

function PlotImage(hObject, eventdata, handles)
% hObject    handle to Save_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Save_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Save_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveAs_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAs_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ScreenResX_Callback(hObject, eventdata, handles)
% hObject    handle to ScreenResX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScreenResX as text
%        str2double(get(hObject,'String')) returns contents of ScreenResX as a double


% --- Executes during object creation, after setting all properties.
function ScreenResX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScreenResX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScreenResY_Callback(hObject, eventdata, handles)
% hObject    handle to ScreenResY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScreenResY as text
%        str2double(get(hObject,'String')) returns contents of ScreenResY as a double


% --- Executes during object creation, after setting all properties.
function ScreenResY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScreenResY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rectROI.
function rectROI_Callback(hObject, eventdata, handles)
% hObject    handle to rectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rectROI


% --- Executes on button press in UseGrid.
function UseGrid_Callback(hObject, eventdata, handles)
% hObject    handle to UseGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseGrid



function RectPos_Callback(hObject, eventdata, handles)
% hObject    handle to RectPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RectPos as text
%        str2double(get(hObject,'String')) returns contents of RectPos as a double


% --- Executes during object creation, after setting all properties.
function RectPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RectPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in polygROI.
function polygROI_Callback(hObject, eventdata, handles)
% hObject    handle to polygROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of polygROI



function roiLabel_Callback(hObject, eventdata, handles)
% hObject    handle to roiLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roiLabel as text
%        str2double(get(hObject,'String')) returns contents of roiLabel as a double


% --- Executes during object creation, after setting all properties.
function roiLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%------------------------------------------------
  
function FigureButtonDownFcn(hObject,eventdata)

ptSelectionTol = .1;
handles = guidata(hObject);

setappdata(hObject,'buttonState',1)

if gco ~= handles.axes1
    return
end


currpt = get(handles.axes1,'currentPoint');
setappdata(hObject,'BtnDownPoint',currpt);

if strcmp(get(hObject,'selectiontype'),'alt')

    setappdata(hObject,'buttonNumber',2)
   
else

    setappdata(hObject,'buttonNumber',1)

end
   
activeroi =  getappdata(hObject,'activeRoi')

roiData = getappdata(hObject,'roiData');

if activeroi ~= 0
    
    xminmax = [min(roiData.instances(activeroi).pos(:,1)) max(roiData.instances(activeroi).pos(:,1)) ]';
    yminmax = [min(roiData.instances(activeroi).pos(:,2)) max(roiData.instances(activeroi).pos(:,2)) ]';
    
    corners = [xminmax([1 2 2 1]), yminmax([1 1 2 2])]; 
    
    insiderect = inpoly(currpt(1),currpt(1),corners(:,1),corners(:,2));
    
    crnorm = diag([diff(xlim) diff(ylim)].^-1); 
    [mindists,activecorner] = min(sqrt(sum( (crnorm*corners- repmat(crnorm*currpt,size(corners,1),1)).^2,2)));
    
    if mindists <= ptSelectionTol
        pivot = corners( mod(acvtivecorner+ 1,4)+1,:);
    else
        pivot = [];
        activecorner=0;
    end
    
    pos = get(handles.axes1,'position');
    
    fignormrect = [corners(1,:) diff(corners(1:2,1)), diff(corners(2:3,2))]
    if acivecorner~=0 

        
        
%         while getappdata(hObject,'buttonState') 
% 
%              updatePlot(hObject,eventdata,roiData.instances(activeroi),pivot,currpt,handles.axes1);
% 
%         end
    end
end    
    
%-----------------------------------------------

function updatePlot(hObject,eventdata,roiInstance,pivot,btnDownPt,axh)
    

    delta = getappdata(axh,currpt)
    roiInstance.pos(:,1:2)

    
%------------------------------------------------
  
function FigureButtonUpFcn(hObject,eventdata)

setappdata(hObject,'buttonState',0)
setappdata(hObject,'buttonNumber',0)
