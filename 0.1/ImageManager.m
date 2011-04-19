function varargout = ImageManager(varargin)
% IMAGEMANAGER M-file for ImageManager.fig
%      IMAGEMANAGER, by itself, creates a new IMAGEMANAGER or raises the existing
%      singleton*.
%
%      H = IMAGEMANAGER returns the handle to a new IMAGEMANAGER or the handle to
%      the existing singleton*.
%
%      IMAGEMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEMANAGER.M with the given input arguments.
%
%      IMAGEMANAGER('Property','Value',...) creates a new IMAGEMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageManager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageManager_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help ImageManager

% Last Modified by GUIDE v2.5 01-Jan-2008 15:03:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageManager_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageManager_OutputFcn, ...
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


% --- Executes just before ImageManager is made visible.
function ImageManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageManager (see VARARGIN)


% Choose default command line output for ImageManager
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
parent = varargin{1};
setappdata(handles.figure1,'parent', parent);
UpdateImage(hObject,eventdata,handles)
setappdata(varargin{1},'activeControl','ImageManager');
setappdata(handles.imageList,'activeImages',[])


setappdata(parent,'ImageManager',handles.figure1);
%     setappdata(handles.figure1,'ImageManagerFunctions',ImageManagerFunctions);
children = getappdata(parent,'children');
children(end+1) = handles.figure1;
setappdata(parent,'children',children)

imageList_Callback(hObject, eventdata, handles);

%subfunctions that can be useful elsewhere
% imFunctions.UpdateImage = @UpdateImage;
% imFunctions.UpdatePosition = @UpdatePosition;
imFunctions.UpdateImage = @(varargin)UpdateImage([],[],handles,varargin{:});
imFunctions.UpdatePosition = @(varargin)UpdatePosition([],[],handles,varargin{:});
imFunctions.LoadImages = @(varargin)loadImages([],[],handles,varargin{:});

setappdata(parent,'ImageManagerFunctions',imFunctions);

MakeActive(hObject, eventdata, handles)

% --- Outputs from this function are returned to the command line.
function varargout = ImageManager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% varargout{2} = imFunctions;

function xpos1_Callback(hObject, eventdata, handles)
% hObject    handle to xpos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xpos1 as text
%        str2double(get(hObject,'String')) returns contents of xpos1 as a double
MakeActive(hObject, eventdata, handles)
UpdatePosition(hObject,eventdata,handles);

UpdateImage(hObject,eventdata,handles);
 activefigure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function xpos1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xpos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xpos2_Callback(hObject, eventdata, handles)
% hObject    handle to xpos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xpos2 as text
%        str2double(get(hObject,'String')) returns contents of xpos2 as a double
MakeActive(hObject, eventdata, handles)
UpdatePosition(hObject,eventdata,handles);

UpdateImage(hObject,eventdata,handles);
 activefigure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function xpos2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xpos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ypos1_Callback(hObject, eventdata, handles)
% hObject    handle to ypos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ypos1 as text
%        str2double(get(hObject,'String')) returns contents of ypos1 as a double
MakeActive(hObject, eventdata, handles)
UpdatePosition(hObject,eventdata,handles);

UpdateImage(hObject,eventdata,handles);
 activefigure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function ypos1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ypos1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ypos2_Callback(hObject, eventdata, handles)
% hObject    handle to ypos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ypos2 as text
%        str2double(get(hObject,'String')) returns contents of ypos2 as a double

MakeActive(hObject, eventdata, handles)
UpdatePosition(hObject,eventdata,handles);

UpdateImage(hObject,eventdata,handles);
 activefigure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function ypos2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ypos2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xres_Callback(hObject, eventdata, handles)
% hObject    handle to xres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MakeActive(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');
screenData = getappdata(parent,'screenData');
imageData = getappdata(parent,'imageData');
CurrentImage = getappdata(parent,'CurrentImage');
CurrentDataSet = getappdata(parent,'CurrentDataSet');

if isempty(CurrentDataSet) || CurrentDataSet == 0, CurrentDataSet = 1; end

newxres = str2double(get(handles.xres,'string'));
if ~isempty(newxres)
    screenData(CurrentDataSet).res(1) = newxres;    
    if CurrentImage~=0
        imageData.images(CurrentImage).screenres = screenData(CurrentDataSet).res;
        setappdata(parent,'imageData',imageData)
    end
    setappdata(parent,'screenData',screenData);
    UpdatePosition(hObject, eventdata, handles)
    UpdateImage(hObject, eventdata, handles)
else
    set(handles.xres,'string',screenData(CurrentDataSet).res(1));
end



% Hints: get(hObject,'String') returns contents of xres as text
%        str2double(get(hObject,'String')) returns contents of xres as a double


% --- Executes during object creation, after setting all properties.
function xres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yres_Callback(hObject, eventdata, handles)
% hObject    handle to yres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');
screenData = getappdata(parent,'screenData');
imageData = getappdata(parent,'imageData');
CurrentImage = getappdata(parent,'CurrentImage');
CurrentDataSet = getappdata(parent,'CurrentDataSet');
if isempty(CurrentDataSet) || CurrentDataSet == 0, CurrentDataSet = 1; end

newyres = str2double(get(handles.yres,'string'));
if ~isempty(newyres)
    screenData(CurrentDataSet).res(2) = newyres;
    
    if CurrentImage~=0
        imageData.images(CurrentImage).screenres = screenData(CurrentDataSet).res;
        setappdata(parent,'imageData',imageData)
    end

    setappdata(parent,'screenData',screenData);
    UpdatePosition(hObject, eventdata, handles)
    UpdateImage(hObject, eventdata, handles)

else
    set(handles.yres,'string',screenData(CurrentDataSet).res(2));
end



% --- Executes during object creation, after setting all properties.
function yres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MakeActive(hObject, eventdata, handles)

% --------------------------------------------------------------------
function LoadImage_Callback(hObject, eventdata, handles)
% Loads a set of images and displays one of them on the main axes
% hObject    handle to LoadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
parent = getappdata(handles.figure1,'parent');

% parent_handles = guidata(parent);

imfmts = imformats;
imexts = strcat('*.',[imfmts.ext]);


% currentTrial = getappdata(parent,'CurrentTrial');
currentImage = getappdata(parent,'CurrentImage');
imageData = getappdata(parent,'imageData');

% dlgans ='Yes';
% if  currentImage>0 && ~isempty(imageData.images(currentImage).filename) 
%     dlgans = questdlg(sprintf('Replace the current image\nfor trial %i with a new one?',currentTrial));
% end
% 
% if ~strcmp(dlgans,'Yes')
%     return;
% end


currpath = '';
if currentImage > 0 && ~isempty(imageData.images(currentImage).filename)
    currpath = imageData.images(currentImage).path;
elseif currentImage > 1 && ~isempty(imageData.images(currentImage-1).filename)
    currpath = imageData.images(currentImage-1).path;
end   
   


[imgs, pth] = uigetfile( {sprintf('%s;',imexts{:}),'Image Files'},'Select an Image File', currpath ,'multiselect','on');

    
if isnumeric(imgs )
    return
end

loadImages(hObject,eventdata,handles,pth,imgs)

%-----------------------------------

function loadImages(hObject,eventdata,handles,pth,imgs)

%Loads a cell array of image files from the specified path. The images
%themselves aren't actually loaded, rather a new image structure array 
% is created.

parent = getappdata(handles.figure1,'parent');


currentImage = getappdata(parent,'CurrentImage');
imageData = getappdata(parent,'imageData');

[imrangenorm, screenres] = GetImagePosition(hObject,eventdata,handles);

imstruct= makeImageStruct(pth, imgs);

for i = 1:length(imstruct)
    imstruct(i).info = imfinfo(strcat(imstruct(i).path,filesep,imstruct(i).filename));
    imstruct(i).xySize = [imstruct(i).info.Width,imstruct(i).info.Height];  
    imstruct(i).code = imageData.codeincr + i;
end    

imageData.codeincr = imstruct(end).code;

irncell = repmat({imrangenorm},1,length(imstruct));
scrncell = repmat({screenres},1,length(imstruct));

multassign = @ (varargin) varargin{:};
[imstruct.position] = multassign(irncell{:});
[imstruct.screenres] = multassign(scrncell{:});

if currentImage == 0 && length(imageData.images)==1 && isempty(imageData.images(1).filename)
    currentImage = 1;
    imageData.images = imstruct;
else
    currentImage = length(imageData.images)+1;

    imageData.images = cat(2,imageData.images,imstruct);
end
    
setappdata(parent,'imageData',imageData)
setappdata(parent,'CurrentImage',currentImage)

UpdateImage(hObject,eventdata,handles)

% set(handles.imageList,'Value',1)
imageList_Callback(hObject, eventdata, handles);
%-----------------------------
function [imrangenorm,screenres] = GetImagePosition(hObject,eventdata,handles)
%Returnsimage Position

MakeActive(hObject, eventdata, handles)

imrange =str2num([get(handles.xpos1,'string'),' ', get(handles.xpos2,'string'),...
        ' ', get(handles.ypos1,'string'), ' ', get(handles.ypos2,'string')]);     %#ok<ST2NM>
screenres =str2num([get(handles.xres,'string'),' ', get(handles.yres,'string')]);
imrangenorm = imrange./screenres([1 1 2 2]);


%-----------------------------
function UpdateImage(hObject,eventdata,handles,varargin)
%
%Redraws the picture on the main axis

if ~ishandle(handles.figure1)
    parent = varargin{1};
else
    parent = getappdata(handles.figure1,'parent');
end

MakeActive(hObject, eventdata, handles) %Make sure the image Manager is active

parent_handles = guidata(parent);
% axes(parent_handles.axes1);
% axis([0 1 0 1]);

currentImage = getappdata(parent,'CurrentImage');
CurrentDataSet = getappdata(parent,'CurrentDataSet');
if isempty(CurrentDataSet) || CurrentDataSet == 0, CurrentDataSet = 1; end



imageData = getappdata(parent,'imageData');
screenData = getappdata(parent,'screenData');

if CurrentDataSet > length(screenData)
    screenData(CurrentDataSet) = screenData(1);
end

screenres = screenData(CurrentDataSet).res;
% screenres =str2num([get(handles.xres,'string'),' ', get(handles.yres,'string')]);
axis(parent_handles.axes1,[0 screenres(1) 0 screenres(2)]);
axis(parent_handles.axes1,'ij');


if currentImage == 0 || currentImage > length(imageData.images) ||isempty(imageData.images(currentImage).filename)
    set(handles.none,'string','none selected','visible','on')
    set(handles.imFileName,'string','','visible','off')
    if ishandle(getappdata(parent,'imageHandle'))
       delete(getappdata(parent,'imageHandle'));
    end
%     activefigure(handles.figure1)
%     axes(parent_handles.axes2);
    set(handles.none,'visible','off');

    return
end


% %Use trial specific screen resolution (although in general this probably
% %won't change, it might might as well be flexible).
% screenData(CurrentDataSet).res = imageData.images(currentImage).screenres;
% screenres = screenData(CurrentDataSet).res;
% 
% setappdata(parent,'screenData',screenData)


% if get(handles.scaleCheckBox,'value')
%     imrange = [0 1 0 1].*screenres([1 1 2 2]);
% elseif get(handles.pixelScaleCheckBox,'value')
%     xysize = imageData.images(currentImage).xySize;
%     screenres = imageData.images(currentImage).screenres;  
%     pos = imageData.images(currentImage).position;    
%     xycenter = ( pos([1 3]) +pos([2 4]))./2.*screenres;
%     imrange = [0 xysize(1) 0 xysize(2)] + xycenter([1 1 2 2]) - .5*xysize([1 1 2 2]);    
% else
%     imrange = imageData.images(currentImage).position.*screenres([1 1 2 2]);
% end

imrange = imageData.images(currentImage).position.*screenres([1 1 2 2]);

im = imread( fullfile(imageData.images(currentImage).path,imageData.images(currentImage).filename));

if size(im,3) == 1
    im = repmat(im,[1,1,3]); %Convert grayscale to truecolor
end

if ishandle(getappdata(parent,'imageHandle'))
    delete(getappdata(parent,'imageHandle'))
end

% imageData.images(currentImage).xySize = [size(im,2) size(im,1)]; 
% imrangenorm = imrange./screenres([1 1 2 2]);

if ishandle(handles.figure1)
    set(handles.xpos1,'String',num2str(imrange(1)))
    set(handles.xpos2,'String',num2str(imrange(2)))
    set(handles.ypos1,'String',num2str(imrange(3)))
    set(handles.ypos2,'String',num2str(imrange(4)))
set(handles.xres,'String',num2str(screenres(1)))
set(handles.yres,'String',num2str(screenres(2)))
set(handles.imFileName,'string',imageData.images(currentImage).filename,'visible','on');
 end

imh = image(imrange(1:2),imrange(3:4), im,'parent',parent_handles.axes1);

axis(parent_handles.axes1,[0 screenres(1) 0 screenres(2)]);

setappdata(parent,'imageHandle',imh);
% axes(parent_handles.axes2);

% activefigure(handles.figure1)




% --- Executes on selection change in imageList.
function imageList_Callback(hObject, eventdata, handles)
% hObject    handle to imageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns imageList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imageList

% pause(.2) %Pause to avoid interrupint a double click
% if ~isequal(get(handles.figure1,'selectiontype'),'open')
%     return
% end

% listItems = get(handles.imageList,'String');
selected = get(handles.imageList,'Value');

if isempty(selected)
    setappdata(handles.imageList,'activeImages',[])
end
% 
% if selected(1) == 1 && length(selected) == 1
%     LoadImage_Callback(hObject, eventdata, handles);
%     return
% elseif selected(1) == 1 
%     selected(1) = [];
% end

% if selected(1) == 2 && length(selected) == 1
if ~isempty(selected) && selected(1) == 1 && selected(1) == 2 
    selected(1) = [];
end


parent = getappdata(handles.figure1,'parent');



imageData = getappdata( getappdata(handles.figure1,'parent'),'imageData');
imgs ={imageData.images.filename};
set(handles.imageList,'String',{'no image',imgs{:}});

if ~isempty(gcbo) && gcbo == handles.imageList && length(selected) == 1 
    setappdata(parent,'CurrentImage',selected(1) - 1)
    UpdateImage(hObject,eventdata,handles);
elseif length(selected) == 1
    set(handles.imageList,'value',getappdata(parent,'CurrentImage')+1);
end

 setappdata(handles.imageList,'activeImages',selected - 1 );
 
 if get(handles.pixelScaleCheckBox,'value')
%      set(handles.pixelScaleCheckBox,'value',0);
     pixelScaleCheckBox_Callback(hObject, eventdata, handles);
 end
 
 if get(handles.scaleCheckBox,'value')
%      set(handles.scaleCheckBox,'value',0);
     scaleCheckBox_Callback(hObject, eventdata, handles);
 end
 activefigure(handles.figure1)
   
% --- Executes during object creation, after setting all properties.
function imageList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function imgListDo_Callback(hObject, eventdata, handles)
% hObject    handle to imgListDo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
% --------------------------------------------------------------------
function RemoveImage_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)
activeImages = getappdata(handles.imageList,'activeImages');

imageData = getappdata(getappdata(handles.figure1,'parent'),'imageData');

if activeImages(1) == 0 && length(activeImages) ==1
    return;
elseif activeImages(1) == 0
   activeImages(1)=[];
end
   

imageData.images( :,activeImages) = [];

parent = getappdata(handles.figure1,'parent');
setappdata(parent ,'imageData',imageData);
setappdata(parent ,'CurrentImage',activeImages(1)-1);

set(handles.imageList,'value', activeImages(1))
imageList_Callback(hObject, eventdata, handles);
 
UpdateImage(hObject,eventdata,handles);
 activefigure(handles.figure1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeActive(hObject, eventdata, handles)
% Makes the current control active in the main window

% if gcf~=handles.figure1
if gcbo~=handles.figure1
    return
end

parent = getappdata(handles.figure1,'parent');

setappdata(parent,'activeControl','ImageManager');

% --- Executes on button press in scaleCheckBox.
function scaleCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to scaleCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.scaleCheckBox,'value')
    set(handles.pixelScaleCheckBox,'value',0)
end

MakeActive(hObject, eventdata, handles)


activeim = getappdata(handles.imageList,'activeImages');

for i = 1:length(activeim)
    UpdatePosition(hObject,eventdata,handles,activeim(i));    % Apply to all selected images
end

UpdateImage(hObject,eventdata,handles);
activefigure(handles.figure1)

 % --- Executes on button press in pixelScaleCheckBox.
function pixelScaleCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to pixelScaleCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


MakeActive(hObject, eventdata, handles)

if get(handles.pixelScaleCheckBox,'value')
    set(handles.scaleCheckBox,'value',0)
end

activeim = getappdata(handles.imageList,'activeImages');

for i = 1:length(activeim)
    UpdatePosition(hObject,eventdata,handles,activeim(i));    % Apply to all selected images
end


UpdateImage(hObject,eventdata,handles);
activefigure(handles.figure1)



function UpdatePosition(hObject, eventdata, handles,currentImage)
% hObject    handle to scaleCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Set image position and scaling

MakeActive(hObject, eventdata, handles)

parent = getappdata(handles.figure1,'parent');

if nargin < 4
    currentImage = getappdata(parent,'CurrentImage');
end

if currentImage == 0
    return
end

imageData = getappdata(parent,'imageData');

if get(handles.scaleCheckBox,'value')
    imrangenorm = [0 1 0 1];
    set(handles.xpos1,'Enable','off')
    set(handles.xpos2,'Enable','off')
    set(handles.ypos1,'Enable','off')
    set(handles.ypos2,'Enable','off')
elseif get(handles.pixelScaleCheckBox,'value')
    xysize = imageData.images(currentImage).xySize;
    screenres = imageData.images(currentImage).screenres;    
    pos = imageData.images(currentImage).position;    
    xycenter = ( pos([1 3]) +pos([2 4]))./2.*screenres;
    imrange = [0 xysize(1) 0 xysize(2)] + xycenter([1 1 2 2]) - .5*xysize([1 1 2 2]);    
    imrangenorm = imrange./screenres([1 1 2 2]);

else
    imrangenorm = imageData.images(currentImage).position;
%     imrangenorm = GetImagePosition(hObject,eventdata,handles);
    set(handles.xpos1,'Enable','on')
    set(handles.xpos2,'Enable','on')
    set(handles.ypos1,'Enable','on')
    set(handles.ypos2,'Enable','on')   
end

imageData.images(currentImage).position = imrangenorm ;
setappdata(parent,'imageData',imageData);


% --------------------------------------------------------------------
function fromText_Callback(hObject, eventdata, handles)
% hObject    handle to fromText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MakeActive(hObject, eventdata, handles)

[file,path] = uigetfile;
if isnumeric(file)
    return
end

imgs = uiImportText(strcat(path,file));
if isempty(imgs)
    return
end


loadImages(hObject,eventdata,handles,path,imgs);



% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fromMatFile_Callback(hObject, eventdata, handles)
% hObject    handle to fromMatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


MakeActive(hObject, eventdata, handles)

[file,path] = uigetfile('*.mat');
if isnumeric(file)
    return
end

imgs = uiImportMatVar(strcat(path,file),'cell');
if isempty(imgs)
    return
end


loadImages(hObject,eventdata,handles,path,imgs);


% --------------------------------------------------------------------
function fromWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to fromWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


MakeActive(hObject, eventdata, handles)


imgs = uiImportWS('cell');

if isempty(imgs)
    return
end

loadImages(hObject,eventdata,handles,cd,imgs);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');

if ishandle(parent)
    setappdata(parent,'ActiveControl','main');
end
