function FX = fixtract(ASL,segment,trialStartTs,trialEndTs,trialIdentifier,paramStruc,Screen,etfs,etparams)

%function FixDat = fixtract(ASL,segment,trialStartT,paramStruc)
%Extracts fixations from ASL data sctructure using velocity criterion

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


%######### PARAMETERS ############
PixelRes = [1024 768];

if nargin < 7 | isempty(Screen)
    screenSize_inches = 17;
    DistancetoScreen_inches = 20; %in inches
else
    DistancetoScreen_inches = Screen.distanceIn; %in inches
    screenSize_inches = Screen.sizeIn;
 
end
if nargin < 8
    etfs = 119.65; %sampling frequency of the eyetracker
end

if nargin < 9
    minFixDur = 100; %Minimum duration of a fixation im ms
    eyeVCriterion = 50 ; %Criterion for saccades in degrees per second.
    RejectFixPtsV = 50 ; %Reject Fixation points if velocity exceeds this value
    pupilRange = [20 36];  %Pupil sizes outside this range are assumed to represent bad discrimination
    ETrangeH = [0 260 ]; %Horizontal screen size in eyetracker coordingates
    ETrangeV = [0 290 ]; %Vertical screen size in eyetracker coordingates
    calibPtSpacing = [408, -304];
    MaxBlinkDur = 200; %Blink duration in milliseconds
     fixptp = .25; %Proportion of time points which must be identified as fixation
else
    minFixDur = etparams.minFixDur; %Minimum duration of a fixation im ms
    eyeVCriterion = etparams.eyeVCriterion; %Criterion for saccades in degrees per second.
    RejectFixPtsV = etparams.RejectFixPtsV; %Reject Fixation points if velocity exceeds this value
    pupilRange = etparams.pupilRange ;  %Pupil sizes outside this range are assumed to represent bad discrimination
    ETrangeH = etparams.ETrangeH ; %Horizontal screen size in eyetracker coordingates
    ETrangeV = etparams.ETrangeV; %Vertical screen size in eyetracker coordingates
    calibPtSpacing =  etparams.calibPtSpacing ;
    MaxBlinkDur = etparams.MaxBlinkDur ; %Blink duration in milliseconds
    fixptp = etparams.fixptp; 
end
    trialPad =200; %Discard fixations starting more than trialPad msec before end of trial
%Define calibration points
screenCenter = PixelRes./2;
calibPts = [kron(calibPtSpacing(1)*[1 1 1]',[-1 0 1]') + screenCenter(1),kron(calibPtSpacing(2)*[-1 0 1]',[1 1 1]')+ screenCenter(2)];
%################################

params = who;  % We will keep all parameters defined above
params(ismember(params,{'ASL','paramStruc'})) = [];

FixDat.aslFileInfo = ASL.fileInfo;


if nargin == 6 & ~isempty(paramStruc)
    for i = 1:length(params)  
        q = paramStruc.(params{i});
        eval(sprintf('%s = q;',params{i}));
    end
end


if nargin < 4 | isempty(trialStartTs) | isempty(trialEndTs)
    trialStartT = 0 ;
    trialEndT = Inf;
end

if nargin < 5 | isempty(trialIdentifier) 
    trialIdentifier = 1:length(trialStartT);
end


for i = 1:length(params)    
    FixDat.parameters.(params{i}) = eval(params{i});
end
%Conversion from screen pixels to angles
inPerPix = screenSize_inches./sqrt(sum(PixelRes.^2));
DegPerPix = atan(inPerPix./DistancetoScreen_inches)*180/pi; %Degrees per pixel


%Get Calibration Data in eyetracker coordinates
for i = 1:size(calibPts,1)

    ETCalibPts(i,1) = str2num(ASL.allHeaderData.Calibration_Values.(sprintf('htgt_data_%i',i)));
    ETCalibPts(i,2) = str2num(ASL.allHeaderData.Calibration_Values.(sprintf('vtgt_data_%i',i)));
end

scrncal = cat(2,calibPts,ones(size(calibPts,1),1));
etcal = cat(2,ETCalibPts,ones(size(ETCalibPts,1),1));
%Transformation matrix from eyetracker to screen pixel coordinates
%ET2SCRN = (calibPts'*etcal)*(etcal'*etcal)^-1; 
ET2SCRN = (scrncal'*etcal)*(etcal'*etcal)^-1; 

minFixPts = round(minFixDur./1000*etfs); %Minimum fixation duration in sampling points

if nargin < 2 & length(ASL.seg) > 1
    segment = input(sprintf('\nThere are %i segments. Which one do you want? ',length(ASL.seg)));
end

etdatH = ASL.seg(segment).horz;
etdatV = ASL.seg(segment).vert;
XDAT = ASL.seg(segment).xdat;

xdatTs = find(diff([0;XDAT]) );
xdatIDs = XDAT(xdatTs);

xdatTs(xdatIDs == 0) = [];
xdatIDs(xdatIDs == 0) = [];

FixDat.xdat.startT = xdatTs'./etfs;
FixDat.xdat.id = xdatIDs'; 


%Time points to discard due to out of range pupil or gaze coordinates
ETbadtimes = ASL.seg(segment).pupil > pupilRange(2) | ASL.seg(segment).pupil < pupilRange(1) ;
ETbadtimes = ETbadtimes | etdatH > ETrangeH(2) | etdatH  < ETrangeH(1);
ETbadtimes = ETbadtimes | etdatV > ETrangeV(2) | etdatV  < ETrangeV(1);

% Extracting short dropout (which may or may not represent blinks)
ShortDrop = zeros(size(ETbadtimes));
BTstart = find(diff([0;ETbadtimes;0]) == 1); 
BTend = find(diff([0;ETbadtimes(1:end-1);0]) == -1); 
badIntvl  = BTend - BTstart; %Duration of bad intervals.
ShortDrop(BTstart(badIntvl <= MaxBlinkDur*etfs/1000)) = 1;
ShortDrop(BTend(badIntvl <= MaxBlinkDur*etfs/1000)) = -1;
ShortDrop = cumsum(ShortDrop);

ETbadtimes(find(ShortDrop)) = 0;
%Eye velocity in deg per second
eyeV = sqrt(sum( ( ET2SCRN(1:2,1:2)^-1 * [0 0;diff(etdatH), diff(etdatV)]').^2 ))'./DegPerPix;

%Smooth velocity
eyeV = convn(eyeV,ones(etparams.smoothwin,1)./etparams.smoothwin,'same');

%Flag saccades wherever velocity exceeds threshold
sacLocs = false(size(eyeV));
fixLocs = false(size(eyeV));
sacLocs(eyeV > eyeVCriterion  & ~ETbadtimes) = true; %time points considered part of saccade if the exceed EyeVCriterion
fixLocs(eyeV < RejectFixPtsV  & ~ETbadtimes & ~ShortDrop) = true; %time points considered part of fixation if they fall below RejectFixPtsV
                                                  %Ambiguous points are
                                                  %rejected at the edges of
                                                  %fixations but included
                                                  %when they occur inside a
                                                  %fixation
                                                  
%Making a structure FixDat which contains information about every fixation.
sacTimes = [find(sacLocs); length(eyeV)];

FixStarts = find(diff(sacTimes) > minFixPts); %Onset of fixation defined by difference of the end of last saccade
                                              %and start of subsequent saccade being greater than threshold  

tr = 0;


FixDat.fix = struct;
count = 1;
disc = [];
for i = 1:length(FixStarts)

         interSacPts = sacTimes(FixStarts(i)):sacTimes(FixStarts(i)+1)-1;
%         timepts = min(interSacPts(fixLocs(interSacPts))):max(interSacPts(fixLocs(interSacPts))); 
         timepts = interSacPts;
         
         if sum(fixLocs(timepts))./(diff(timepts([1 end]))+1) < fixptp;
             disc(end+1) = FixStarts(i);
             continue;
             
         end
        tr = find(sacTimes(FixStarts(i)) > trialStartTs & sacTimes(FixStarts(i)) < trialEndTs - round(trialPad*etfs/1000));
        if ~isempty(tr) & ~isempty(timepts)
            meanPos = [mean(etdatH(timepts)),mean(etdatV(timepts)),1]*ET2SCRN';

            FixDat.fix(count).timepts = uint32(timepts);  %Timepoints contained in fixation
            FixDat.fix(end).meanPos = meanPos(1:2);          %Mean of gaze coordinate during fixation
            

            FixDat.fix(end).var = [var(etdatH(timepts)),var(etdatV(timepts))]; %Variance of gaze coordinates
            if length(FixDat.fix) > 1
                FixDat.fix(end).sacVec = FixDat.fix(end).meanPos - FixDat.fix(end-1).meanPos; %Direction  of saccade at onset of fixation
            else                                                            %relative to previous fixation (not defined for the first one)
                FixDat.fix(end).sacVec = [nan nan]; 
            end

            FixDat.fix(end).startT = timepts(1)./etfs;
%             FixDat.fix(end).startT = timepts(1);
            
            FixDat.fix(end).dur = diff(timepts([1,end]))./etfs;
%             FixDat.fix(end).durT = diff(timepts([1,end]));
            %Fixation is suspect if it starts during or immediately after pupil
            %loss or abberant position data 
            FixDat.fix(end).suspect = sum(ETbadtimes(interSacPts(1)+[-1:0])) > 0;   
            % A dispersion of  more than 2 degrees is also suspect      
            if sqrt(sum((ET2SCRN(1:2,1:2)*sqrt(FixDat.fix(end).var)').^2)) > 2./DegPerPix
                FixDat.fix(end).suspect = 1;
            end     
            xdindex = max(find(xdatTs < timepts(1)));
            currxdat = xdatIDs(xdindex);
            if isempty(currxdat)
                currxdat = 0;
                xdindex = 0;
            end
                
            if xdindex > 1
                FixDat.fix(end).xdhist = [xdatIDs( xdindex - 1), currxdat];
            
            else
                FixDat.fix(end).xdhist = [0, currxdat];
            end
            FixDat.fix(end).xdat = currxdat;
            FixDat.fix(end).xdindex  = xdindex ;
            
            FixDat.fix(end).trial = trialIdentifier(tr);
            count = count+1;
        end
        
end

FixDat.etfs= etfs;

FX.et2scrn = ET2SCRN;
FX.units = 'xy Screen';
FX.seg = FixDat;
FX.seg = FixDat;

%FixDat.fix(find([FixDat.fix.suspect])) = [];


