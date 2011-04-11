function imdata = loadImages(h,imagefiles,fpath)


%  imdata = loadImages(h,imagefiles,fpath)
%
%  Loads specified image files.
%
%
% See also ImageManager



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end


if nargin < 2
    imagefiles = [];
end
if nargin < 3
    fpath = [];
end


handles = guidata(h);

GazeReader('Image_Manager_Callback',h,[],handles);
imgfun  = getappdata(h,'ImageManagerFunctions');
imgfun.LoadImages(fpath,imagefiles);

if nargout > 0
    imdata = getappdata(h,'imageData');
end

