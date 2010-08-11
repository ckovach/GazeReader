function importData(h,datafiles,fpath)

%  importData(h,datafiles)
%
% Imports eye-tracker data files. Currently supported file types are Eyelink (*.edf) and ASL (*.asl).
%

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

GazeReader('importEytrackerData_menu_Callback',h,[],handles,datafiles,fpath);