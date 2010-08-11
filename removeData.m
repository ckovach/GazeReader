function removeData(h,datindices)

%  removeData(h,dataindices)
%
% Removes data from the currently loaded stack.
%
% See also importMatData and importData



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end

ehd = getappdata(h,'eyetrackerHeaderData');

if nargin < 2
    dataindices = 1:length(ehd);
end

handles = guidata(h);

dsfun = getappdata(h,'DataSetWindowFunctions');

dsfun.delete(dataindices);


