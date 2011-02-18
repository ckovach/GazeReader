function varargout = basisFcnDlg(varargin)
% BASISFCNDLG M-file for basisFcnDlg.fig
%      BASISFCNDLG, by itself, creates a new BASISFCNDLG or raises the existing
%      singleton*.
%
%      H = BASISFCNDLG returns the handle to a new BASISFCNDLG or the handle to
%      the existing singleton*.
%
%      BASISFCNDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BASISFCNDLG.M with the given input arguments.
%
%      BASISFCNDLG('Property','Value',...) creates a new BASISFCNDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before basisFcnDlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to basisFcnDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help basisFcnDlg

% Last Modified by GUIDE v2.5 18-Feb-2011 10:56:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @basisFcnDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @basisFcnDlg_OutputFcn, ...
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


% --- Executes just before basisFcnDlg is made visible.
function basisFcnDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to basisFcnDlg (see VARARGIN)

parent = varargin{1};
setappdata(hObject,'parent',parent);

set(handles.basis_menu,'String',varargin{2})


setappdata(parent,'basisSet',[]);
setappdata(parent,'basisOrd',[]);

% Choose default command line output for basisFcnDlg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes basisFcnDlg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = basisFcnDlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in basis_menu.
function basis_menu_Callback(hObject, eventdata, handles)
% hObject    handle to basis_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns basis_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from basis_menu


% --- Executes during object creation, after setting all properties.
function basis_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basis_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function order_num_Callback(hObject, eventdata, handles)
% hObject    handle to order_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of order_num as text
%        str2double(get(hObject,'String')) returns contents of order_num as a double


% --- Executes during object creation, after setting all properties.
function order_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');

setappdata(parent,'basisSet',get(handles.basis_menu,'value'));
setappdata(parent,'basisOrd',str2num(get(handles.order_num,'string')));
setappdata(parent,'basisKeepDC',get(handles.keepdc,'value'));

delete(handles.figure1)




% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% pushbutton1_Callback(hObject, eventdata, handles)
% cancel_Callback(hObject, eventdata, handles)


% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent = getappdata(handles.figure1,'parent');


setappdata(parent,'basisSet',[]);
setappdata(parent,'basisOrd',[]);

delete(handles.figure1)


% --- Executes on button press in keepdc.
function keepdc_Callback(hObject, eventdata, handles)
% hObject    handle to keepdc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keepdc
