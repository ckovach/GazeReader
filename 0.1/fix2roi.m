
function [FixRoi] = fix2roi(FixDat,stimdat,varargin)

% Maps fixation onto regions of interest
%
% Input:
% 
%    stimdat - A 1xN structure, where N is the number of trials. The structure
%           must contain at least 1 field which defines the position of 
%           each of M ROIs. This can take the form of a Mx2 vector with xy
%           or ij coordinates (default is xy). Alternately, the field can
%           be a 1xM structure with fields giving the vertices of
%           polygonal ROIs. In the former case, a second field must contain
%           either a scalar defining a radius from the position parameter
%           of points within a circle around the position or a vector
%           defining the vertices of the ROI polynomial.
%
% FixDat  - Fixation position starting time and duration.
%
%Output:
%   FixRoi  - Fixation with respect to regions of interest and the classification of the Roi.

% RoiUnits = 'ij'; %Is roi position expressed in xy or ij units?

% roi_center_offset = [0 0];
trialOnsets = zeros(size(stimdat));
trialnumfield = 'xdat';
durfield = 'dur';
fxposfield = 'meanPos';
shiftfield = 'shiftvec';
roiradius = 110;

% if isfield(FixDat,'parameters')
%     ScreenParam = FixDat.parameters;
% end


% fcpos = [];
useshapefield = 1;
i = 1;
xdattrialval = [];
trials  = [];
while i <= length(varargin)
   switch lower(varargin{i})
       case 'roiposfield'
          roiposfield = varargin{i+1};
            i = i+1;
       case 'roishapefield'
          roishapefield = varargin{i+1};
            i = i+1;
       case 'roiradius'
          roiradius = varargin{i+1};
            i = i+1;
            useshapefield = 0;
%        case 'roicenteroffset'
%             roi_center_offset  = varargin{i+1};
%             i = i+1;
       case 'trialonsets'
               trialOnsets  = varargin{i+1};
            i = i+1;
%        case 'roiunits'
%            RoiUnits = varargin{i+1};
%             i = i+1;
       case 'trialnumfield'
           trialnumfield = varargin{i+1};
            i = i+1;
       case 'trials' %Trials data supplied in a separate vector
            trials = varargin{i+1};
            i = i+1;
           
%        case 'screenparam'
%            ScreenParam = varargin{i+1};
%             i = i+1;
       case 'durfield'
           durfield = varargin{i+1};
            i = i+1;
       case 'fxposfield'
           fxposfield = varargin{i+1};
            i = i+1;
       case 'shiftfield'
           shiftfield = varargin{i+1};
            i = i+1;
       case 'xdathigh'  %Xdat remains high for the duration of the trial with one of the given values  
           xdattrialval = varargin{i+1};
            i = i+1;

%        case 'fixcrosspos'
%           fcpos = varargin{i+1};
%             i = i+1;
        otherwise

           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end


% if isempty(fcpos)
%     FixCrossPos = ScreenParam.screenCenter;
% else
%    FixCrossPos = fcpos;
% end

if isempty(xdattrialval) && isempty(trials)
    trials = [FixDat.fix.(trialnumfield)];
elseif isempty(trials)
    trialid = [FixDat.fix.(trialnumfield)];
    trialon = ismember(trialid ,xdattrialval);
    trials = trialon .* cumsum(diff([0,trialid]) ~= -[0,trialid(1:end-1)] & diff([0,trialid]) ~= 0 );
end    

%for i = 1:length(trialOnsets)
ntrials = length(stimdat);


% if strcmp(RoiUnits,'ij');
%     unitCorrOrder = [2 1];
%     unitCorrOffset = ScreenParam.PixelRes(2);
%    
% else
%     unitCorrOrder = [1 2];
% end

%function handles for simple utilities
cellmeanfun = @(x) mean(x,1);
getcol1 = @(M) M(:,1);
getcol2 = @(M) M(:,2);



for i = 1:ntrials
%     
%     if i == 18
%         keyboard
%     end
    getfix = find(trials == i);
    if ~isempty(getfix)
        fix = cat(1,FixDat.fix(getfix).(fxposfield))'; %Position of fixations in this trial
        
        if isstruct(stimdat(i).(roiposfield)) 
            %Distance to ROI is computed wih respect to the CENTER- this
            %may not be appropriate for highly irregular polygons!
            cm = cellfun(cellmeanfun,{stimdat(i).(roiposfield).(roishapefield)});
            
            Rois = cat(1,cm{:})';   %Position of ROIs
        else
            Rois = stimdat(i).(roiposfield)';
        end
        
        
        %Matrix of distances to each ROI for each fixation
        Ds = permute(sqrt(sum((Rois(:,:,ones(1,size(fix,2))) - permute(fix(:,:,ones(1,length(Rois))),[1 3 2])).^2,1)),[2 3 1]);
        
        %Matrix of angles to each ROI for each fixation in clockwise
        %radians with 0 deg being up.
        AnglesToRoi = permute( atan2( Rois(1,:,ones(1,size(fix,2))) - permute(fix(1,:,ones(1,length(Rois))),[1 3 2]) , Rois(2,:,ones(1,size(fix,2))) - permute(fix(2,:,ones(1,length(Rois))),[1 3 2])  )  ,[2 3 1]);
        
        [mn,mi] = min(Ds);
        mindists = zeros(size(Ds));
        mindists(mi + (0:size(Ds,2)-1)*size(Ds,1)) = 1; 
        
        if ~isstruct(stimdat(i).(roiposfield)) 
           
            % Fixation is associated with ROI if it falls within the specified
            % radius and is the minimum.
            if useshapefield && length(stimdat(i).(roishapefield)) == 1 
                [targ,q] = find(Ds < stimdat(i).(roishapefield) & mindists );                
            elseif useshapefield 
                
                %Modify for fixed polygonal
                %[targ,q] = find(Ds < stimdat(i).(roishapefield) & mindists );                
                error('This part of the function is under construction. To be continued...')
            else
                [targ,q] = find(Ds < roiradius & mindists );                
            end                
                
        else
            inpolyfun = @(A,B) inpolygon(fix(:,1),fix(:,2),A,B) ;
            
            Cx = cellfun(getcol1,{stimdat(i).(roiposfield).(roishapefield)},'UniformOutput',false);
            Cy = cellfun(getcol2,{stimdat(i).(roiposfield).(roishapefield)},'UniformOutput',false);
            
            [targ,q] = find(cellfun(inpolyfun, Cx,Cy)  &  mindists );        
        end
        
        FixRoi(i).roi = zeros(1,size(fix,2));
        FixRoi(i).em = zeros(1,size(fix,2));
        FixRoi(i).id = zeros(1,size(fix,2));
        FixRoi(i).roi(q) = targ;            %Sequence of targets
        FixRoi(i).roiShift = (diff([0,FixRoi(i).roi]) ~=0);  %1 for shifts between ROIs, 0 otherwise
        FixRoi(i).roiOverlap =  sum(Ds < roiradius) > 1; %1 if fixation is lands in overlapping rois. Only the roi with 
                                                             % the minimum distance is reported in the roi field

        % When ROIs overlap, report alternate rois in the altRoi field  FIX
        % FOR POLYGONAL ROIS
       FixRoi(i).altRoi = cell(1,length(FixRoi(i).roi));
       for ro = find(FixRoi(i).roiOverlap)
            FixRoi(i).altRoi{ro} =  find(Ds(:,ro) < roiradius);  
            %FixRoi(i).altRoi{ro}(FixRoi(i).roi(ro) == FixRoi(i).altRoi{ro}) = [];
       end
        
        %Compute number of fixations on a given target
        RoiBlock = FixRoi(i).roi;
        RoiBlock(~FixRoi(i).roiShift) = -1; %Count only fixations that shift from one Roi to another
        RoiBlockSum = cumsum(RoiBlock(ones(1,size(Rois,2)+1),:) == (0:size(Rois,2))'*ones(1,size(fix,2)),2);
        FixRoi(i).fixRep = RoiBlockSum(FixRoi(i).roi+1 + (0:length(FixRoi(i).roi)-1)*size(RoiBlockSum,1)); %NUmber of times the current ROI has been fixed, not counting repeats
        FixRoi(i).fixNum = cumsum(RoiBlock >= 0); %Overall order of fixation, not including repeats
        FixRoi(i).fixRoiNum = cumsum(RoiBlock > 0);%Order of fixation excluding fixations at ROI 0 (outside any specified ROIs).  
        FixRoi(i).PointNum = 1:size(fix,2);%Order of data point in trial.  
        FixRoi(i).startT = [FixDat.fix(getfix).startT]; %Fixation onset times
        FixRoi(i).trialT = [FixDat.fix(getfix).startT] - trialOnsets(i); %Onset with respect to trial onset
         FixRoi(i).trial = i;
        FixRoi(i).(durfield) = [FixDat.fix(getfix).(durfield)]';  %Fixation duration
        FixRoi(i).fixdat = FixDat.fix(getfix);  %Original fixation data associated with this trial
        FixRoi(i).fixIndex = getfix;  %Indices in the FIX structure
        
        
%         if FixRoi(i).fixdat(1).dt - FixRoi(i).trialT(1) < .750  %Assume fixation is at fixation cross when dt exceeds 750 msec        
             initfix = FixRoi(i).fixdat(1).(fxposfield) - FixRoi(i).fixdat(1).(shiftfield); %fixation point before the onset of the trial
%         else
%             initfix = FixCrossPos;
%         end
        
        InitD = sqrt(sum((Rois - initfix( ones(size(Rois,2),1),:)').^2));  %Initial distances to current fixation at trial onset
        [s,smat] = sort(cat(2,InitD' ,Ds));  %Sorting distances

        InitAng = atan2( Rois(1,:) - initfix( ones(size(Rois,2),1),1)' , Rois(2,:) - initfix( ones(size(Rois,2),1),2)');
%         rsmat = smat;

        %Rearranging ranks into the order of ROIs (so that row i corresponds to ROI
        %i again)
        rsmat = zeros(size(smat)); 
        rsmat(smat + size(smat,1)*ones(size(smat,1),1)*(0:size(smat,2)-1)) = (1:size(smat,1))'*ones(1,size(smat,2));
        rsmat(:,find(FixRoi(i).roi)+1) = rsmat(:,find(FixRoi(i).roi)+1)  - 1; %Subtract 1 from Colums which contain a roi

 %         Dmat = sqrt(sum(diff([ FixCrossPos;cat(1,FixDat.fix(getfix).(fxposfield))]',1,2).^2));
        Dmat = sqrt(sum(diff([ initfix;cat(1,FixDat.fix(getfix).(fxposfield))],1,1)'.^2)); %Distance between fixations
        FixRoi(i).Dmat = Dmat; 
        FixRoi(i).Drank = -1*ones(size(Dmat));
        FixRoi(i).Drank(FixRoi(i).roi ~= 0) = rsmat(FixRoi(i).roi(FixRoi(i).roi ~=0 ) + (find(FixRoi(i).roi) - 1)*size(Rois,2));
        FixRoi(i).DmatAllRois = cat(2,InitD',Ds(:,1:end-1)); %Distance to all ROIs
        FixRoi(i).AngleAllRois = cat(2, InitAng', AnglesToRoi(:,1:end-1)); 
        FixRoi(i).DrankAllRois = rsmat(:,1:end-1); 
         FixRoi(i).startpos = initfix; %Posiion of the last fixation before trial onset if interfixation time is less than 750 msec,
                                        %otherwise, the location of the fixation cross.
    else
         FixRoi(i).trial = i;
    end
end

