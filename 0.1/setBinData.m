function setBinData(h,bindata)


%  setBinData(h,bindata)
%
% Associates the bindata with current data set
% and brings up the bin data window
%
% bindata is a structure returned by MAKEBINDATA
%
% See also makeBinData

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end



if nargin > 1
    setappdata(h,'binData',bindata);
end

handles = guidata(h);

GazeReader('binManagerMenu_Callback',h,[],handles);

