function model = setModelData(h,modeldata)


%  setRegData(h,modeldata)
%
% Associates the modeldata with current data set
% and brings up the modelManager window
%
% modeldata.models(i) are structures returned by makeModel
%



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end



if nargin > 1
    setappdata(h,'modelData',modeldata);
end

handles = guidata(h);

GazeReader('models_menu_Callback',h,[],handles);

% if ~isempty(tmfun)
%     tmfun.update();
% end
if nargout > 0
    modeldata = getappdata(h,'modelData');
end






