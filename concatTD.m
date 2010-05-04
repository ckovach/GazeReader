function [TD,lasttrs] = concatTD(lastendts,lastfxns,lastBGRs,xdmaps,varargin)

%TDcat = concatTD(TD1,TD2,TD3,...)
%Concatenates multiple trialData structures into a single one.


if length(varargin)>1
    TDin = cat(2,varargin{:});
else
    TDin = varargin{1};
end


lastendt = 0;
lastcode = 0;
lastnumber = 0;
lastfxn = 0;
lastBGR = 0;
lasttrn = 0;

TD = TDin(1);
TD.trials= cat(2,TDin.trials);
for i = 1:length(TDin)
    if ~isempty(TDin(i).trials)
%         if i == 1
%             TD= TDin(i);
%         else

            
             for k = 1:length(TDin(i).trials)   

%                  TD.trials(end + k) = TDin(i).trials(k);
                 TD.trials(lasttrn + k).number = TDin(i).trials(k).number + lastnumber;
                 TD.trials(lasttrn + k).startTime = TDin(i).trials(k).startTime + lastendt;
                 TD.trials(lasttrn + k).stopTime = TDin(i).trials(k).stopTime + lastendt;
                 TD.trials(lasttrn + k).fixOnsetTimes = TDin(i).trials(k).fixOnsetTimes + lastendt;
                 TD.trials(lasttrn + k).fixations = TDin(i).trials(k).fixations + lastfxn;
%                  TD.trials(lasttrn + k).fixOnsetTimes = TDin(i).trials(k).fixOnsetTimes + lastendt;
                 TD.trials(lasttrn + k).samplePts = TDin(i).trials(k).samplePts + lastendt;
                 TD.trials(lasttrn + k).code = TDin(i).trials(k).code + lastcode;
                 TD.trials(lasttrn + k).binGroup = TDin(i).trials(k).binGroup + lastBGR;
                 TD.trials(lasttrn + k).startCode = xdmaps{i}(ismember(xdmaps{i},TDin(i).trials(k).startCode));
                 TD.trials(lasttrn + k).stopCode = xdmaps{i}(ismember(xdmaps{i},TDin(i).trials(k).stopCode));
             end
           
          
%         end
        
        lastendt = lastendts(i);
        lastBGR = lastBGRs(i);
        lastcode = TD.trials(lasttrn + k).code;
        lastnumber = TD.trials(lasttrn + k).number;
        lastfxn = lastfxns(i);
        lasttrn = lasttrn + length(TDin(i).trials);
        lasttrs(i) = lastnumber;
        TD.codeincr = lastcode;
    end
end

        