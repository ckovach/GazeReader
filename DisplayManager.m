function varargout = DisplayManager(varargin)
% DISPLAYMANAGER M-file for DisplayManager.fig
%      DISPLAYMANAGER, by itself, creates a new DISPLAYMANAGER or raises the existing
%      singleton*.
%
%      H = DISPLAYMANAGER returns the handle to a new DISPLAYMANAGER or the handle to
%      the existing singleton*.
%
%      DISPLAYMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAYMANAGER.M with the given input arguments.
%
%      DISPLAYMANAGER('Property','Value',...) creates a new DISPLAYMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DisplayManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DisplayManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DisplayManager

% Last Modified by GUIDE v2.5 07-Jan-2008 22:17:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DisplayManager_OpeningFcn, ...
                   'gui_OutputFcn',  @DisplayManager_OutputFcn, ...
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


% --- Executes just before DisplayManager is made visible.
function DisplayManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DisplayManager (see VARARGIN)

% Choose default command line output for DisplayManager
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DisplayManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DisplayManager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in fixationList.
function fixationList_Callback(hObject, eventdata, handles)
% hObject    handle to fixationList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fixationList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fixationList


% --- Executes during object creation, after setting all properties.
function fixationList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixationList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fixEdit_Callback(hObject, eventdata, handles)
% hObject    handle to fixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixEdit as text
%        str2double(get(hObject,'String')) returns contents of fixEdit as a double


% --- Executes during object creation, after setting all properties.
function fixEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in jumpFixUp.
function jumpFixUp_Callback(hObject, eventdata, handles)
% hObject    handle to jumpFixUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in jumpTrialUp.
function jumpTrialUp_Callback(hObject, eventdata, handles)
% hObject    handle to jumpTrialUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in jumpTrialDown.
function jumpTrialDown_Callback(hObject, eventdata, handles)
% hObject    handle to jumpTrialDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in jumpFixDown.
function jumpFixDown_Callback(hObject, eventdata, handles)
% hObject    handle to jumpFixDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


