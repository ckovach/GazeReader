function [FIX,RAW] = readEDF(filename,nhist)


%Uses edfmex to read edf files and convert them to a format consistent with
%that used for ASL files. 

if nargin < 2 
    nhist = 2;
end

if ~exist(filename,'file')
    error('Unable to find: %s',filename);
end

try 
    EDF = edfmex(filename);
catch
    error(sprintf(['\nIt appears the pre-compiled file edfmex.mex is broken or incompatible witht your system.',...
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

code = {EDF.FEVENT.codestring};
msg = strcmp(code,'MESSAGEEVENT');
msgevnt = EDF.FEVENT(msg);
msgcodes = {msgevnt.message};
msgt = [msgevnt.sttime] - EDF.FSAMPLE.time(1)+1;

[unqmsg,q,msgmap] = unique(msgcodes);

msgsig  = zeros(size(EDF.FSAMPLE.time));

msgsig(msgt) = sum(repmat(msgt,length(msgt),1) == repmat(msgt',1,length(msgt))) ; 
xdatraw = msgmap(cumsum(msgsig));

if nargout > 1
    RAW.segData.horz = double(EDF.FSAMPLE.gx(eye,:));
    RAW.segData.vert = double(EDF.FSAMPLE.gy(eye,:));
    RAW.segData.pupil = double(EDF.FSAMPLE.pa(eye,:));
    RAW.segData.xdat = xdatraw;
    RAW.segData.degConversion(:,1) = double(EDF.FSAMPLE.rx);
    RAW.segData.degConversion(:,2) = double(EDF.FSAMPLE.ry);
    RAW.segData.time = double(EDF.FSAMPLE.time);
    RAW.segData.xdatCodes = unqmsg;
end

%fixation data
endfix = strcmp(code,'ENDFIX');
endsacc = strcmp(code,'ENDSACC');


fxs = EDF.FEVENT(endfix);
if length(msgevnt) ~= length(msgcodes)
    warning(sprintf('Uh Oh. Some of the messages fields are empty.\nEverything will be out of alignment if I continue.'))        
    inp = input(sprintf('1. Stop\n2. Continue\n: '));
    if inp == 1
         error('Stopping...');
    end

end

for i = 1:length(msgmap)
    FIX.seg.xdat(i).startT = msgevnt(i).sttime;
    FIX.seg.xdat(i).id = msgmap(i);
    FIX.seg.xdat(i).code= msgcodes{i};    
end
msgts = [FIX.seg.xdat.startT];

for i = 1:length(fxs)
    FIX.seg.fix(i).meanPos = [fxs(i).gavx fxs(i).gavy];
    FIX.seg.fix(i).startT = fxs(i).sttime;
    FIX.seg.fix(i).dur   = fxs(i).entime -fxs(i).sttime; 
    FIX.seg.fix(i).endT = fxs(i).entime;
    FIX.seg.fix(i).updconv = ([fxs(i).supd_x, fxs(i).supd_y] + [fxs(i).eupd_x, fxs(i).eupd_y])./2;
    FIX.seg.fix(i).eye = fxs(i).eye;

    recentmsg = msgts <= fxs(i).sttime;
    frecentmsg= find(recentmsg );
    FIX.seg.fix(i).xdhist  = zeros(1,nhist);
    FIX.seg.fix(i).xdhist(end:-1:end-min([length(frecentmsg) nhist]-1)) = msgmap(frecentmsg(end:-1:end-min([length(frecentmsg) nhist]-1)));
    FIX.seg.fix(i).xdat = FIX.seg.fix(i).xdhist(end);
    cr = cumsum(recentmsg);
    FIX.seg.fix(i).xdindex = cr(end);    
    FIX.seg.fix(i).lastcode = msgcodes{frecentmsg(end)};
    

    if i>1
        FIX.seg.fix(i).shiftvec = diff(cat(1,FIX.seg.fix(i + [-1:0]).meanPos));
        FIX.seg.fix(i).dt = diff(cat(1,FIX.seg.fix(i + [-1:0]).startT)); % difference between fixation onsets of 
                                                                                         %the current and previous fixation
    else
        FIX.seg.fix(i).shiftvec =[nan nan];
        FIX.seg.fix(i).dt =nan;
    end
      FIX.units = 'xy Screen';        
end

saccs = EDF.FEVENT(endsacc);

for i = 1:length(saccs)
    FIX.seg.sac(i).startT= saccs(i).sttime;
    FIX.seg.sac(i).dur= saccs(i).entime-saccs(i).sttime;
    FIX.seg.sac(i).startPos= [saccs(i).gstx saccs(i).gsty];
    FIX.seg.sac(i).shiftVec= [saccs(i).gstx saccs(i).gsty] - [saccs(i).genx saccs(i).geny];    
    
    recentmsg = msgts <= saccs(i).sttime;
    frecentmsg = find(msgts <= saccs(i).sttime);
    FIX.seg.sac(i).xdhist  = zeros(1,nhist);
    FIX.seg.sac(i).xdhist(end:-1:end-min([length(frecentmsg) nhist]-1)) = msgmap(frecentmsg(end:-1:end-min([length(frecentmsg) nhist]-1)));
    FIX.seg.sac(i).xdindex = cumsum(recentmsg);
    FIX.seg.sac(i).lastcode = msgcodes{frecentmsg(end)};
    
end
FIX.seg.xdatCodes = unqmsg;

