function importMatData(h,datafiles,fpath)

%  importMatData(h,datafiles)
%
% Imports a mat file containing FIX and RAW structures returned by ReadEDF
%
% See also READEDF

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end

if nargin < 2
    datafile = [];
elseif ~iscell(datafiles)
    datafiles = {datafiles};
end

if nargin < 3
    fpath = [];
end

handles = guidata(h);

GazeReader('importMatFile_menu_Callback',h,[],handles,datafiles,fpath);