
function matvarout = uiImportText(matfile,classtype)

% lines = importText(txtfile)
%   A simple utility to select a variable from a mat file
%   and import it based on the users selection.



    
    fig = figure;

    h = uicontrol('Parent',fig,'Style','listbox','units','normalized','position',[.1 .1 .8 .8]);
    uicontrol('Parent',fig,'Style','text','String','Select lines...','fontsize',12,'units','normalized','position',[.1 .9 .3 .05]);



    rg = {vars.name};
    size = {vars.size};

    set(h,'String',rg,'Callback',@Callback,'max',1);

    setappdata(h,'lines',rg)
    % iskb = 0;
    set(fig,'currentcharacter',char(1))

    
%     get(fig,'selectiontype');
%     while ( isempty(double(get(fig,'currentcharacter'))) || double(get(fig,'currentcharacter')) ~= 13) && ~strcmp(get(fig,'selectiontype'),'open')
    while  ~strcmp(get(fig,'selectiontype'),'open')

        uiwait

    end

    selectedLines = getappdata(h,'selectedlines');

    delete(fig)

 ld = load(matfile,selectedLines);
 matvarout = ld.(selectedLines);
 
%%%%%%%%%%%%%%%%%%%

function Callback(hObject,eventData)

selected = get(hObject,'value');

rg = getappdata(hObject,'lines');

setappdata(hObject,'selectedlines',rg(selected));

uiresume

