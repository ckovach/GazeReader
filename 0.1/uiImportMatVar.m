
function matvarout = uiImportMatVar(matfile,classtype)

% lines = importText(txtfile)
%   A simple utility to select a variable from a mat file
%   and import it based on the users selection.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------



vars = whos('-file',matfile);
if nargin > 2
    getclass = strcmp({vars.class},classtype);
    vars(~getclass) = [];
else
    classtype  ='';
end

matvarout= [];

if isempty(vars)    
    warning('No %s variables found in file!',classtype)
    return
end

if length(vars)>1
    
    fig = figure;
    set(fig,'menubar','none');

    h = uicontrol('Parent',fig,'Style','listbox','units','normalized','position',[.1 .1 .8 .8]);
    uicontrol('Parent',fig,'Style','text','String','Select variable...','fontsize',12,'units','normalized','position',[.1 .9 .3 .05]);



    rg = {vars.name};
    size = {vars.size};

    set(h,'String',rg,'Callback',@Callback,'max',1);

    setappdata(h,'lines',rg)
    % iskb = 0;
    set(fig,'currentcharacter',char(1))
   
   set(fig,'selectiontype','normal');

    
%     get(fig,'selectiontype');
    while ( isempty(double(get(fig,'currentcharacter'))) || double(get(fig,'currentcharacter')) ~= 13) && ~strcmp(get(fig,'selectiontype'),'open')
        uiwait
        if ~ishandle(fig)
              return
        end
    end

    selectedLines = getappdata(h,'selectedlines');

    delete(fig)
    
else
   
    selectedLines = {vars.name};
end
ld = load(matfile,selectedLines{1});
matvarout = ld.(selectedLines{1});
 
%%%%%%%%%%%%%%%%%%%

function Callback(hObject,eventData)

selected = get(hObject,'value');

rg = getappdata(hObject,'lines');

setappdata(hObject,'selectedlines',rg(selected));

uiresume

