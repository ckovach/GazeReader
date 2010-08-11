function regdata = setRegData(h,regdata)


%  setRegData(h,regdata)
%
% Associates the regdata with current data set
% and brings up the reg data window
%
% regdata.regressors(i) are structures returned by makeregressor
%



if nargin < 1 || isempty(h) || ~ishandle(h)    
    h = GazeReader;
end



if nargin > 1
    setappdata(h,'regData',regdata);
end

handles = guidata(h);

GazeReader('Regressors_menu_Callback',h,[],handles);

% if ~isempty(tmfun)
%     tmfun.update();
% end
if nargout > 0
    regdata = getappdata(h,'regData');
end






