function [FIX,RAW] = readEDF(filename,varargin)

% [FIX,RAW] = readEDF(filename,varargin)
%
%Uses edfmex to read edf files and convert them to a format consistent with
%that used for ASL files. 
%
% See also EDFMEX

nhist = 2;
i = 1;
subtractFirstXdat = false;
eyeinp = 'x';
while i <= length(varargin)
   switch lower(varargin{i})
       case 'subtract_first_xdat_time' %Subtract the time of the first message event. 
           subtractFirstXdat = true;
       case 'eye' %'Use EYE: 'L' for left,'R' for Right
           eyeinp = varargin{i+1};
           i = i+1;
       case 'nhist' %'Use EYE: 'L' for left,'R' for Right
           eyeinp = varargin{i+1};
           i = i+1;
        otherwise
           error([varargin{i},' is not a valid option.']);
   end         
    i = i+1;
end

% if nargin < 2 || isempty(nhist)
%     nhist = 2;
% end

if ~exist(filename,'file')
    error('Unable to find: %s',filename);
end

try 
    fprintf('%s ',filename)
    fprintf('\n')
    EDF = edfmex(filename);
    fprintf('\n')
catch
    le = lasterror;
    error(sprintf(['\nIt appears the pre-compiled file edfmex.mex is broken or incompatible with your system.',...
            '\nTry recompiling edfmex.cpp in the source directory with mex ''edfmex.cpp edfapi.lib''']))
end

RAW.recdata = EDF.RECORDINGS;
RAW.eyetracker = 'Eyelink';
RAW.header = EDF.HEADER;
RAW.filename = EDF.FILENAME;
FIX.recdata = EDF.RECORDINGS;
FIX.eyetracker = 'Eyelink';
FIX.header = EDF.HEADER;
FIX.filename = EDF.FILENAME;

eye = EDF.RECORDINGS(1).eye;
etfs = EDF.RECORDINGS(1).sample_rate;
if eye == 3
%     inp = 'x';
    while ~ismember(eyeinp,{'R','L'})
        eyeinp = inputdlg('Which eye to analyze: L or R?');
    end
    
    eye = isequal(eyeinp{1},'L')+1;
end
    
code = {EDF.FEVENT.codestring};
msg = strcmp(code,'MESSAGEEVENT');
msgevnt = EDF.FEVENT(msg);
msgcodes = {msgevnt.message};
msgt = [msgevnt.sttime] - EDF.FSAMPLE.time(1)+1;

[unqmsg,q,msgmap] = unique(msgcodes);

msgsig  = zeros(size(EDF.FSAMPLE.time));

msgsig(msgt) = sum(repmat(msgt,length(msgt),1) == repmat(msgt',1,length(msgt))) ; 
xdatraw = msgmap(cumsum(msgsig));

if subtractFirstXdat
    startTime = msgevnt(1).sttime;
else
    startTime = 0;
end

if nargout > 1
    RAW.seg.horz = double(EDF.FSAMPLE.gx(eye,:));
    RAW.seg.vert = double(EDF.FSAMPLE.gy(eye,:));
    RAW.seg.pupil = double(EDF.FSAMPLE.pa(eye,:));
    RAW.seg.xdat = xdatraw;
    RAW.seg.degConversion(:,1) = double(EDF.FSAMPLE.rx);
    RAW.seg.degConversion(:,2) = double(EDF.FSAMPLE.ry);
    RAW.seg.time = double(EDF.FSAMPLE.time) - double(startTime) ;
    RAW.seg.xdatCodes = 1:length(unqmsg);
%     RAW.seg.xdatLabels = unqmsg;
end

%fixation data
endfix = strcmp(code,'ENDFIX');
endsacc = strcmp(code,'ENDSACC');


fxs = EDF.FEVENT(endfix);

if length(msgevnt) ~= length(msgcodes)
    warning('Uh Oh. Some of the messages fields are empty.\nEverything will be out of alignment if I continue.')        
    inp = input(sprintf('1. Stop\n2. Continue\n: '));
    if inp == 1
         error('Stopping...');
    end

end

for i = 1:length(msgmap)
    FIX.seg.xdat(i).startT = double(msgevnt(i).sttime) - double(startTime);
    FIX.seg.xdat(i).id = msgmap(i);
    FIX.seg.xdat(i).code= msgcodes{i};    
end
msgts = [FIX.seg.xdat.startT];

saccs = EDF.FEVENT(endsacc);
if ~isempty(saccs)
    saceye = double([saccs.eye]);
    saccs = saccs(saceye == eye -1);
    sact = [saccs.sttime];
end

    
i = 0;
for k = 1:length(fxs)
    if double(fxs(k).eye) ~= eye -1
        continue
    end
    i = i+1;
    
    FIX.seg.fix(i).meanPos = [fxs(k).gavx fxs(k).gavy];
    FIX.seg.fix(i).startT = double(fxs(k).sttime-startTime);
    FIX.seg.fix(i).dur   = fxs(k).entime -fxs(k).sttime; 
    FIX.seg.fix(i).endT = double(fxs(k).entime-startTime);
    FIX.seg.fix(i).updconv = ([fxs(k).supd_x, fxs(k).supd_y] + [fxs(k).eupd_x, fxs(k).eupd_y])./2;
    FIX.seg.fix(i).eye = fxs(k).eye;

    recentmsg = msgts <= fxs(k).sttime-startTime;
    frecentmsg= find(recentmsg );
    FIX.seg.fix(i).xdhist  = zeros(1,nhist);
    FIX.seg.fix(i).xdhist(end:-1:end-min([length(frecentmsg) nhist]-1)) = msgmap(frecentmsg(end:-1:end-min([length(frecentmsg) nhist]-1)));
    FIX.seg.fix(i).xdat = FIX.seg.fix(i).xdhist(end);
    cr = cumsum(recentmsg);
    FIX.seg.fix(i).xdindex = cr(end);    
    FIX.seg.fix(i).lastcode = msgcodes{frecentmsg(end)};
    
    if i>1
        FIX.seg.fix(i).shiftvec = diff(cat(1,FIX.seg.fix(i + (-1:0)).meanPos));
        FIX.seg.fix(i).dt = diff(cat(1,FIX.seg.fix(i + (-1:0)).startT)); % difference between fixation onsets of 
    
        FIX.seg.fix(i).sac = find( sact <= FIX.seg.fix(i).startT,1,'last');
        
                                                                                         %the current and previous fixation
    else
        FIX.seg.fix(i).shiftvec =[nan nan];
        FIX.seg.fix(i).dt =nan;
    end
      FIX.units = 'xy Screen';        
end


i = 0;
fxts = [FIX.seg.fix.startT];
for k = 1:length(saccs)
    if double(saccs(k).eye) ~= eye-1
        continue
    end
    i = i+1;
    
    stt = double(saccs(k).sttime);
    dfst = fxts - stt;
    lastfix = find(dfst < 0,1,'last'); 
    nextfix = find(dfst > 0,1,'first'); 
    
    FIX.seg.sac(i).startT= stt;
    FIX.seg.sac(i).dur= saccs(k).entime-saccs(k).sttime;
    FIX.seg.sac(i).startPos= [saccs(k).gstx saccs(k).gsty];
    FIX.seg.sac(i).shiftVec= [saccs(k).gstx saccs(k).gsty] - [saccs(k).genx saccs(k).geny];    
    FIX.seg.sac(i).eye= saccs(k).eye;
    if ~isempty(lastfix)
        FIX.seg.sac(i).lastfix= lastfix;
    else
        FIX.seg.sac(i).lastfix= nan;        
    end
    if ~isempty(nextfix)
        FIX.seg.sac(i).nextfix= nextfix;
    else
        FIX.seg.sac(i).nextfix= nan;        
    end
    recentmsg = msgts <= saccs(k).sttime;
    frecentmsg = find(msgts <= saccs(k).sttime);
    FIX.seg.sac(i).xdhist  = zeros(1,nhist);
    FIX.seg.sac(i).xdhist(end:-1:end-min([length(frecentmsg) nhist]-1)) = msgmap(frecentmsg(end:-1:end-min([length(frecentmsg) nhist]-1)));
    cr = cumsum(recentmsg);
    FIX.seg.sac(i).xdindex = cr(end);
    FIX.seg.sac(i).lastcode = msgcodes{frecentmsg(end)};
    
end
FIX.seg.xdatCodes = unqmsg;
RAW.fs = etfs;
FIX.fs = etfs;
