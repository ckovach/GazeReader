function poolData(h,dataindices)

%  removeData(h,dataindices)
%
% Concatenates specified data sets
%
% See also REMOVEDATA

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end

ehd = getappdata(h,'eyetrackerHeaderData');

if nargin < 2
    dataindices = 1:length(ehd);
end


dsfun = getappdata(h,'DataSetWindowFunctions');

dsfun.pool(dataindices);


