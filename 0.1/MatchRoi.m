
function [RoiDat,FixRoi] = matchRoi(FixDat,imdat,varargin)

% Matches fixation data in FixDat with ROI data in imdat
%
% Summary of data structures used in this script
% 
% imdat  - Stimulus data saved in searchImdata_*
% FixDat - Fixation position starting time and duration.
% RoiDat - Region of interest location and stimulus number with respect to data in imdat.
% FixRoi - Fixation with respect to regions of interest and the classification of the Roi.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

RoiUnits = 'ij'; %Is roi position expressed in xy or ij units?

roi_center_offset = [0 0];
trialOnsets = zeros(size(imdat.trialdat));
trialnumfield = 'trial';
durfield = 'dur';
posfield = 'meanPos';
shiftfield = 'shiftvec';
startTfield = 'startT';ntrials = length(imdat.trialdat);

durationfield = 'dt';
conversion = 1;
if isfield(FixDat,'parameters')
    ScreenParam = FixDat.parameters;
end

ntrials = length(imdat.trialdat);

gettrials = 1:ntrials;

fcpos = [];

i = 1;
while i <= length(varargin)
   switch lower(varargin{i})
       case 'roicenteroffset'
            roi_center_offset  = varargin{i+1};
            i = i+1;
       case 'trialonsets'
               trialOnsets  = varargin{i+1};
            i = i+1;
       case 'roiunits'
           RoiUnits = varargin{i+1};
            i = i+1;
       case 'trialnumfield'
           trialnumfield = varargin{i+1};
            i = i+1;
       case 'screenparam'
           ScreenParam = varargin{i+1};
            i = i+1;
       case 'durfield'
           durfield = varargin{i+1};
            i = i+1;
       case 'posfield'
           posfield = varargin{i+1};
            i = i+1;
       case 'shiftfield'
           shiftfield = varargin{i+1};
            i = i+1;
       case 'fixcrosspos'
          fcpos = varargin{i+1};
            i = i+1;
       case 'conversionfactor'
          conversion = varargin{i+1};
            i = i+1;
       case 'starttfield'
          startTfield = varargin{i+1};
            i = i+1;
       case 'durationfield'
          durationfield = varargin{i+1};
            i = i+1;
       case 'trialnumbers'
          trialnumbers = varargin{i+1};
            i = i+1;
        otherwise

           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end


if isempty(fcpos)
    FixCrossPos = ScreenParam.screenCenter;
else
   FixCrossPos = fcpos;
end

trials = [FixDat.fix.(trialnumfield)];

[emolabel,q,emonum] = unique(imdat.targdat.emo);
[idlabel,q,idnum] = unique(imdat.targdat.id);
RoiRadius = 110;
% RoiDat 
%for i = 1:length(trialOnsets)


if strcmp(RoiUnits,'ij');
    unitCorrOrder = [2 1];
    unitCorrOffset = ScreenParam.PixelRes(2);
   
else
    unitCorrOrder = [1 2];
end

for i = 1:max(trialnumbers)
   %i = trials(trnum);
   if sum(trialnumbers == i) >0
       RoiDat(i).Center = cat(1,imdat.trialdat(trialnumbers == i).pos{:}) + roi_center_offset(ones(1,length([imdat.trialdat(trialnumbers == i).pos])),:);

        if strcmp(RoiUnits,'ij');
            RoiDat(i).Center = RoiDat(i).Center(:,unitCorrOrder);

            RoiDat(i).Center(:,2) =  unitCorrOffset - RoiDat(i).Center(:,2); %Flip vertical axis if ROIs are in ij units
         end

       stim = cat(1,imdat.trialdat(trialnumbers == i).num{:});
       RoiDat(i).stimNum = stim; 
       RoiDat(i).ID = idnum(stim); 
       RoiDat(i).EM = emonum(stim);
       RoiDat(i).rad = RoiRadius ;
   end
   RoiDat(i).trial = i;
end


    
%for i = 1:length(trialOnsets)
for i = 1:max(trialnumbers)
%     
%     if i == 18
%         keyboard
%     end
    getfix = find(trials == i);
    if ~isempty(getfix)
        fix = cat(1,FixDat.fix(getfix).(posfield))'*conversion; %Position of fixations in this trial
        Rois = RoiDat(i).Center';   %Position of ROIs
      
        %Rois = RoiDat(i).Center(:,unitCorrOrder)';   %Position of ROIs
      
        %Matrix of distances to each ROI for each fixation
        Ds = permute(sqrt(sum((Rois(:,:,ones(1,size(fix,2))) - permute(fix(:,:,ones(1,length(Rois))),[1 3 2])).^2,1)),[2 3 1]);
        
        %Matrix of angles to each ROI for each fixation in clockwise
        %radians with 0 deg being up.
        AnglesToRoi = permute( atan2( Rois(1,:,ones(1,size(fix,2))) - permute(fix(1,:,ones(1,length(Rois))),[1 3 2]) , Rois(2,:,ones(1,size(fix,2))) - permute(fix(2,:,ones(1,length(Rois))),[1 3 2])  )  ,[2 3 1]);
        
        % Fixation is associated with ROI if it falls within the specified
        % radius and is the minimum.
        [mn,mi] = min(Ds);
        mindists = zeros(size(Ds));
        mindists(mi + [0:size(Ds,2)-1]*size(Ds,1)) = 1; 
        [targ,q] = find(Ds < RoiDat(i).rad & mindists );

        
        FixRoi(i).roi = zeros(1,size(fix,2));
        FixRoi(i).em = zeros(1,size(fix,2));
        FixRoi(i).id = zeros(1,size(fix,2));
        FixRoi(i).roi(q) = targ;            %Sequence of targets
        FixRoi(i).em(q) = RoiDat(i).EM(targ);   %Expression of target
        FixRoi(i).id(q) = RoiDat(i).ID(targ);   % Identity of target
        FixRoi(i).stimNum(q) = RoiDat(i).stimNum(targ); %Stimulus number
        FixRoi(i).roiShift = (diff([0,FixRoi(i).roi]) ~=0);  %1 for shifts between ROIs, 0 otherwise
        FixRoi(i).roiOverlap =  sum(Ds < RoiDat(i).rad) > 1; %1 if fixation is lands in overlapping rois. Only the roi with 
                                                             % the minimum distance is reported in the roi field

        % When ROIs overlap, report alternate rois in the altRoi field
       FixRoi(i).altRoi = cell(1,length(FixRoi(i).roi));
       for ro = find(FixRoi(i).roiOverlap)
            FixRoi(i).altRoi{ro} =  find(Ds(:,ro) < RoiDat(i).rad);  
            %FixRoi(i).altRoi{ro}(FixRoi(i).roi(ro) == FixRoi(i).altRoi{ro}) = [];
       end
        
        %Compute number of fixations on a given target
        RoiBlock = FixRoi(i).roi;
        RoiBlock(~FixRoi(i).roiShift) = -1; %Count only fixations that shift from one Roi to another
        RoiBlockSum = cumsum(RoiBlock(ones(1,size(Rois,2)+1),:) == [0:size(Rois,2)]'*ones(1,size(fix,2)),2);
        FixRoi(i).fixRep = RoiBlockSum(FixRoi(i).roi+1 + [0:length(FixRoi(i).roi)-1]*size(RoiBlockSum,1)); %NUmber of times the current ROI has been fixed, not counting repeats
        FixRoi(i).fixNum = cumsum(RoiBlock >= 0); %Overall order of fixation, not including repeats
        FixRoi(i).fixRoiNum = cumsum(RoiBlock > 0);%Order of fixation excluding fixations at ROI 0 (outside any specified ROIs).  
        FixRoi(i).PointNum = 1:size(fix,2);%Order of data point in trial.  
        FixRoi(i).startT = [FixDat.fix(getfix).(startTfield)]; %Fixation onset times
        FixRoi(i).trialT = [FixDat.fix(getfix).(startTfield)] - trialOnsets(trialnumbers == i); %Onset with respect to trial onset
         FixRoi(i).trial =i;
        FixRoi(i).(durfield) = [FixDat.fix(getfix).(durfield)]';  %Fixation duration
        FixRoi(i).fixdat = FixDat.fix(getfix);  %Original fixation data associated with this trial
        FixRoi(i).fixIndex = getfix;  %Indices in the FIX structure
        
        CentD = sqrt(sum((Rois - FixCrossPos( ones(size(Rois,2),1),:)').^2));
%         [s,smat] = sort(cat(2,CentD',Ds)); % assumes initial fixation is
%                                           at fixation cross
        
        if FixRoi(i).fixdat(1).(durationfield) - FixRoi(i).trialT(1) < .750  %Assume fixation is at fixation cross when dt exceeds 750 msec        
            initfix = FixRoi(i).fixdat(1).(posfield)*conversion - FixRoi(i).fixdat(1).(shiftfield)*conversion; %fixation point before the onset of the trial
        else
            initfix = FixCrossPos;
        end
        
        InitD = sqrt(sum((Rois - initfix( ones(size(Rois,2),1),:)').^2));  %Initial distances to current fixation at trial onset
        [s,smat] = sort(cat(2,InitD' ,Ds));  %Sorting distances

        InitAng = atan2( Rois(1,:) - initfix( ones(size(Rois,2),1),1)' , Rois(2,:) - initfix( ones(size(Rois,2),1),2)');
%         rsmat = smat;

        %Rearranging ranks into the order of ROIs (so that row i corresponds to ROI
        %i again)
        rsmat = zeros(size(smat)); 
        rsmat(smat + size(smat,1)*ones(size(smat,1),1)*[0:size(smat,2)-1]) = [1:size(smat,1)]'*ones(1,size(smat,2));
        rsmat(:,find(FixRoi(i).roi)+1) = rsmat(:,find(FixRoi(i).roi)+1)  - 1; %Subtract 1 from Colums which contain a roi

 %         Dmat = sqrt(sum(diff([ FixCrossPos;cat(1,FixDat.fix(getfix).(posfield))]',1,2).^2));
        Dmat = sqrt(sum(diff([ initfix;cat(1,FixDat.fix(getfix).(posfield))*conversion]',1,2).^2)); %Distance between fixations
        FixRoi(i).Dmat = Dmat; 
        FixRoi(i).Drank = -1*ones(size(Dmat));
        FixRoi(i).Drank(find(FixRoi(i).roi)) = rsmat(FixRoi(i).roi(find(FixRoi(i).roi)) + (find(FixRoi(i).roi) - 1)*size(Rois,2));
        FixRoi(i).DmatAllRois = cat(2,InitD',Ds(:,1:end-1)); %Distance to all ROIs
        FixRoi(i).AngleAllRois = cat(2, InitAng', AnglesToRoi(:,1:end-1)); 
        FixRoi(i).DrankAllRois = rsmat(:,1:end-1); 
         FixRoi(i).startpos = initfix; %Position of the last fixation before trial onset if interfixation time is less than 750 msec,
                                        %otherwise, the location of the fixation cross.
    else
         FixRoi(i).trial = i;
    end
end

% 0
