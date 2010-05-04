function [FX,trend] =  CorrectFixationDrift(FX, crossintvlmarks, fclocation, priorD, timeD,porder ,segs)

% function [ trend, P] =  CorrectFixationDrift(FX, crossintvlmarks, fcpos, priorCenter, priorD, timeD);
%   FX - fixations 
%   crossintvlevents - xdat markers delimiting fixation cross
%   fcpos - fixation cross position (remember to make this the same units 
%       as FX!)
%   priorD - standard squared width of gaussian weighting function (points farther from the expected 
%       fixation cross location are weighted less, which is a way of discounting abberant and stray fixations) 
%   timeD - An alternative way of specifying when the fixation cross
%       happens: a two element cell array with a vector containing xdat
%       signals
%       that follow FC and a fixed time interval as second element indicating the delay to fc onset.
% Fits a  polynomial to mean fixation position during the fixation cross
% epoch. Allows correction of calibration drift.
%

if nargin < 6 || isempty(porder)
    porder = 4; %Order of drift polynomial
end

if nargin < 7 || isempty(segs)
    segs = 1:length(FX.seg);
end


for s = segs
    
    if size(fclocation,1) == 1

        fcpos = ones(length(FX.seg(s).fix),1)*fclocation;

    else
        fcpos = fclocation;
    end



    fxts = [FX.seg(s).fix.startT];
    fxpos = cat(1,FX.seg(s).fix.meanPos);

    if nargin < 4

        weights = ones(length(fxts),1);
    else

        weights = exp( - diag((fxpos - fcpos)*priorD^-1*(fxpos - fcpos)'));    

    end

%     weights = weights./sum(weights);



    if nargin >= 5
        trialonsetT = [FX.seg(s).xdat(ismember([FX.seg(s).xdat.id], timeD{1})).startT];
        crossonsetT = trialonsetT - timeD{2};
    elseif size(crossintvlmarks,1) == 1
        crossonsetT = [FX.seg(s).xdat([FX.seg(s).xdat.id] == crossintvlmarks(1)).startT];
        trialonsetT = [FX.seg(s).xdat([FX.seg(s).xdat.id] == crossintvlmarks(2)).startT];
    else
        crossonsetT  = crossintvlmarks(:,1);
        trialonsetT = crossintvlmarks(:,2);
    end    

    %fctrialnum = [];
    FCpos = [];
    FCts = [];
    allinds = [];
    for i = 1:length(crossonsetT)

        inds = find(fxts > crossonsetT(i)  & fxts < trialonsetT(i));
        allinds = cat(2,allinds,inds);
        %warning('off','MATLAB:divideByZero')
        %fcavpos(i,:) = mean(fxpos(inds,:),1);
        if ~isempty(inds) 
            fclastfix(i,:) = fxpos(inds(end),:);
        else
            fclastfix(i,:) = [nan nan];
        end
        FCpos = cat(1,FCpos,fxpos(inds,:) - fcpos(inds,:));
        FCts = cat(2,FCts,fxts(inds));    

    end
    mxfc = max(FCts);
    FCts = FCts./max(FCts);
    weights = weights(allinds);
    weights = weights./sum(weights);
    X = FCts(ones(1,porder+1),:).^((ones(length(FCts),1)*(porder:-1:0))');

    P = (X*diag(weights)*FCpos)'*(X*diag(weights)*X')'^-1; %Weighted Fitting of polynomial

    % P(1,:) = polyfit(FCts',FCpos(:,1) ,porder);
    % P(2,:) = polyfit(FCts',FCpos(:,2) ,porder);

    FX.seg(s).trendP = P;


    for i = 1:length(FX.seg(s).fix)
        trend(i,:) = [polyval(P(1,:),FX.seg(s).fix(i).startT./mxfc), polyval(P(2,:),FX.seg(s).fix(i).startT./mxfc)];
        FX.seg(s).fix(i).posdetrend = FX.seg(s).fix(i).meanPos - trend(i,:);
        %FX.seg.fix(i).shiftdetrend = FX.seg.fix(i).shiftvec - diff(trend([i,i-1],:));

    end

end