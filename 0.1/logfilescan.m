function S = logfilescan(logfilename)

%Scans a Presentation (Neurobehavioral systems) log file and returns a
%structure containing the times and data associated with events.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C. Kovach 2007

S = struct('Scenario','','Time',[],'Picture',struct([]),'Response',struct([]),'Sound',struct([]));

raw = textread(logfilename,'%s','delimiter','\n');

decell = @(x)[x{:}];
S.Scenario = decell(regexp(raw{1},'Scenario - (.*)','tokens','once'));
S.Time= decell(regexp(raw{2},'Logfile written - (.*)','tokens','once'));

for i = 1:length(raw)


%        q = regexp(raw{i},'([\S]*)\t*(\d*)\t*(Picture|Response|Sound)\t*(\w[^\t]*\w)\t*(\d*)\t*(\d*)\t*(\d*)\t*(\d*)\t*(\d*)\t*(\d*)[^\n]','tokens')
   q = regexp(raw{i},'([\S]*)\t*(\d*)\t*(Picture|Response|Sound)\t*([^\t]*)\t*(\d*)\t*(\d*)\t*(\d*)\t*(\d*)\t*[^\n]*','tokens');
 
   if ~isempty(q) & ~isempty(q{1}{2})
      
       q = q{1}; 
       S.(q{3})(end+1).code = q{4};
       S.(q{3})(end).trial = str2num(q{2});
       S.(q{3})(end).subject = q{1};
       S.(q{3})(end).time = str2num(q{5});
       S.(q{3})(end).ttime = str2num(q{6});
       S.(q{3})(end).terr = str2num(q{7});       
%        S.(q{3})(end).derr= str2num(q{8});
%        S.(q{3})(end).reqt= str2num(str2num(q{9}));
%        S.(q{3})(end).reqd= str2num(q{10});
%        S.(q{3})(end).reqd= str2num(q{11});
   end
end

