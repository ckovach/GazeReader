

function trialData = makeTrialData(varargin)


% function trialData = makeTrialData('property',value)
% Intitialize a trial sructure.
%   
%           trialData.trials - a structure array with the following fiels:
% 
%         number        -   Trial number
%         startTime     -   start time (ms)
%         stopTime      -   stop time (ms)
%         startCode     -   code for xdat which marks the onset of the trial
%         startCodeLabel-   code label for event which marks the offset of the trial
%         startEventCode-   code for event which marks the onset of the trial
%         stopCode      -   code for xdat which marks the offset of the trial
%         stopCodeLabel -   etc
%         stopEventCode -   etc
%         samplePts     -   samples within the raw data contained in trial
%         info          -
%         image         -   code for the image associated with the trial
%         binGroup      -   vector of codes for bin group(s) associated with the trial
%         roiIndex      -   
%         seg           -   
%         fixations     -   indices of fixations within trial
%         fixmat        -   matrix of fixation x binmembership indicator
%                           variable
%         fixOnsetTimes -   fixation onset times for fixations within the trial
%         fixbin        -   bin numbers which contain the fixations
%         nfix          -   number of fixations within the trial
%         nbin          -   number of bins within the trial
%         binarea       -   areas of bins within the trials
%
% See also SETTRIALDATA


newtrials  = struct('number',[],'startTime',[],'stopTime',[],'startCode',[],'startCodeLabel',[],'startEventCode',[],...
                    'stopCode',[],'stopCodeLabel',[],'stopEventCode',[],'samplePts',[],'info',[],'image',[],'binGroup',[],...
                    'roiIndex',[],'seg',[],'fixations',[],'fixmat',[],'fixOnsetTimes',[],'fixbin',[],'nfix',[],'nbin',[],...
                    'binareas',[],'code',[]);



if nargin > 0 && isstruct(varargin{1})
    oldTrialData = varargin{1};
    varargin(1) = [];
else
    oldTrialData.trials = newtrials([]);
    oldTrialData.codeincr = 0;    
end

if nargin == 0 || isempty(varargin{1})
    trialData = oldTrialData;
    return
end

i = 1;

m2c = @(x) mat2cell(x,1,ones(size(x)));
while i <= length(varargin)
   switch (varargin{i})
       case 'startCode'
%           starrtcode = varargin{i+1};
          i = i+1;
       case 'stopCode'
%           stopcode = varargin{i+1}; %1Xntrials logical vector with 1  whenever the to-be-assigned bin occurs
          i = i+1;
       case 'startTime'
%           starrtime = varargin{i+1}; %1Xntrials logical vector with 1  whenever the to-be-assigned bin occurs
          i = i+1;
       case 'stopTime'
%           stoptime = varargin{i+1}; %Total number of trials
          i = i+1;
       case 'image'
%           imageIndex = varargin{i+1}; %Total number of trials
          i = i+1;
       case 'binGroup'
%           bingrpIndex = varargin{i+1}; %Total number of trials
          i = i+1;
       case 'code'
          i = i+1;
       case 'fixbin'
          i = i+1;
       case 'binGroup'
          i = i+1;
       case 'roiIndex'
%           roiIndex = varargin{i+1}; %Total number of trials
          i = i+1;
       case 'startCodeLabel'
           i = i+1;
       case'stopCodeLabel'
           i = i+1;
       case'startEventCode'
           i = i+1;
       case'stopEventCode'
           i = i+1;
       case'number'
           i = i+1;
       case'nfix'
           i = i+1;
       case'nbin'
           i = i+1;
        otherwise
           error([varargin{i},' is not a valid option.']);
   end         
   
   if isempty(varargin{i})
           c = {[]};       
   elseif ~iscell(varargin{i})
           c = m2c(varargin{i});       
   else
       c = varargin{i};
   end
   
   newtrials(length(c)).(varargin{i-1}) = [];   
   [newtrials.(varargin{i-1})] = c{:};
   
   i = i+1;
   
   
end


 c = m2c(1:length(newtrials));
[newtrials.number] = c{:};


if isempty([newtrials.code])
    c = m2c((1:length(newtrials)) + oldTrialData.codeincr);
    [newtrials.code] = c{:};
end
trialData.trials= cat(2,oldTrialData.trials,newtrials);
trialData.codeincr = max([newtrials.code]);


