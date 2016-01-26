
function EventData = makeEventData(varargin)

% EventData = makeEventData([fieldname],value,...)
% EventData = makeEventData(EventData, [fieldname],value,...)
% 
% Intitializes an EventData Structure, or appends events to an existing
% one, with the following fields:
%
%   .xdat  - structure of inputs recorded by the eye tracker with fields
%              .startT - onset time in ms
%              .id     - id code of event
%              .code   - numerical or character code describing event 
%
%   .event - structure of events with fields
%       .label - label for event
%       .time  - time in ms for event
%       .type  - type of event: 0 - nothing, 1 - trial onset, 2 - trial offset
%       .xdatcide - associated xdatcode
%       .info - miscellaneous info.
%       .code - unique identifying code

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------



EventData.events = struct('label',[],'time',[],'type',[],'xdatcode',[],'info',[],'code',[]);
EventData.events(1) = [];
EventData.xdat = struct('startT',[],'id',[],'code',[]);
EventData.xdat(1) = [];
EventData.codeincr = 0;
codeincr = 0;
events = EventData.events;

if nargin == 0 || isempty(varargin{1})
    return
end

if isstruct(varargin{1})
    EventData = varargin{1};
    if isfield(EventData,'codeincr')
        codeincr = EventData.codeincr;
    end  
    if ~isfield(EventData,'events')
        EventData.events = events;
    end
       
    i = 2;
else
    i = 1;
end
m2c = @(x) mat2cell(x,1,ones(size(x)));

while i <= length(varargin)
   switch (varargin{i})
       case 'label'
          i = i+1;
       case 'time'
          i = i+1;
       case 'type'
          i = i+1;
       case 'code'
          i = i+1;
       case 'info'
          i = i+1;
       case 'xdatcode'
          i = i+1;
        otherwise
           error([varargin{i},' is not a valid option.']);
   end         
   
   if ~iscell(varargin{i})
       c = m2c(varargin{i});
   else
       c = varargin{i};
   end
   events(length(varargin{i})).(varargin{i-1}) = 0;   
   [events.(varargin{i-1})] = c{:};
   
   i = i+1;
   
   
end
if isempty([events.code])     
    c = m2c((1:length(events)) + codeincr);   
   [events.code] = c{:};   
end

EventData.events = cat(2,EventData.events,events);

evts = [EventData.events.time]; 
[srtt,srtind] = sort(evts);

EventData.codeincr = max([EventData.events.code]);

%Sort into temporal orders
EventData.events = EventData.events(srtind);




