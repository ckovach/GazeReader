function makeActiveDataSet(h,dataset)


%  makeActiveDataSet(h,dataset)
%
% Makes the data set with the specified index active

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
    setappdata(h,'CurrentDataSet',dataset);
end


