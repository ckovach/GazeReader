
function selectedLines = uiImportText(txtfile)

% lines = importText(txtfile)
%   A simple utility to display the contents of a text file in a list box
%   and to import the user selected lines.


fig = figure;
set(fig,'menubar','none');

h = uicontrol('Parent',fig,'Style','listbox','units','normalized','position',[.1 .1 .8 .8]);
uicontrol('Parent',fig,'Style','text','String','Select lines...','fontsize',12,'units','normalized','position',[.1 .9 .3 .05]);

fid = fopen(txtfile,'r');
rg = regexp(char(fread(fid)'),'[^\n]*','match');
fclose(fid);
selectedLines = [];

set(h,'String',rg,'Callback',@Callback,'max',1e6);

setappdata(h,'lines',rg)
% iskb = 0;
set(fig,'currentcharacter',char(1))
while ( isempty(double(get(fig,'currentcharacter'))) || double(get(fig,'currentcharacter')) ~= 13)...
        && ~strcmp(get(fig,'selectiontype'),'open')
    
    uiwait
     if ~ishandle(fig)
        return
    end
end

selectedLines = getappdata(h,'selectedlines');

delete(fig)
%%%%%%%%%%%%%%%%%%%

function Callback(hObject,eventData)

selected = get(hObject,'value');

rg = getappdata(hObject,'lines');

setappdata(hObject,'selectedlines',rg(selected));

uiresume

