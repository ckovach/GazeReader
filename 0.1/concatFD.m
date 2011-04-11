function [FD, lastendts, lastfxns,xdmaps] = concatFD(varargin)

%FDcat = concatFD(FD1,FD2,FD3,...)
%Concatenates multiple fixationData structures into a single one.


if length(varargin)>1
    FDin = cat(2,varargin{:});
else
    FDin = varargin{1};
end

tpad = 1e3;

lastendt = 0;
fixns = nan(1,length(FDin));
sacns = nan(1,length(FDin));
xdns = zeros(1,length(FDin));
for i = 1:length(FDin)
    fixns(i) = length(FDin(i).fix);
    if isfield(FDin(i),'sac');
        sacns(i) = length(FDin(i).sac);
    end
   xdns(i) = length(FDin(i).xdat);
end

lastfxn = 0;
lastsacn = 0;
lastxdn = 0;

FD = FDin(1);
FD.fix = cat(2,FDin.fix);
if isfield(FDin,'sac')
    FD.sac = cat(2,FDin.sac);
end

%%FIX XDAT CODE ASSIGNMENT
FD.xdatCodes = unique(cat(2,FDin.xdatCodes));

for i = 1:length(FDin)

%         if i == 1
%             FD= FDin(i);
%         else

            [ism,xdatmap]  = ismember(FDin(i).xdatCodes,FD.xdatCodes);
             for k = 1:length(FDin(i).fix)   

%                  FD.fix(lastfxn + k) = FDin(i).fix(k);
                 FD.fix(lastfxn + k).startT = FDin(i).fix(k).startT + lastendt;
                 FD.fix(lastfxn + k).endT = FDin(i).fix(k).endT + lastendt;
                 if FDin(i).fix(k).xdat ~=0
                     if isnumeric(FDin(i).xdatCodes)    
                         FD.fix(lastfxn + k).xdat = xdatmap(FDin(i).fix(k).xdat == FDin(i).xdatCodes);
                     else
                          FD.fix(lastfxn + k).xdat =xdatmap(FDin(i).fix(k).xdat );
                     end
                 else
                     FD.fix(lastfxn + k).xdat = 0;
                 end
                 xdh = FDin(i).fix(k).xdhist;
                 for q = 1:length(xdh)
                     if xdh(q)~=0
                          if isnumeric(FDin(i).xdatCodes)    
                            FD.fix(lastfxn + k).xdhist(q) = xdatmap(xdh(q)==FDin(i).xdatCodes);
                          else
                            FD.fix(lastfxn + k).xdhist(q) = xdatmap(xdh(q));
                          end                              
                     else
                         FD.fix(lastfxn + k).xdhist(q) = 0;
                     end
                 end
                 FD.fix(lastfxn + k).xdindex = FDin(i).fix(k).xdindex + lastxdn;
                 
             end
          
           for k = 1:length(FDin(i).xdat)   
%                  FD.xdat(lastfxn +  k) = FDin(i).xdat(k);
                 FD.xdat(lastxdn + k).startT = FDin(i).xdat(k).startT + lastendt;
                 if isnumeric(FDin(i).xdatCodes)             
                     FD.xdat(lastxdn + k).id = xdatmap(FDin(i).xdat(k).id == FDin(i).xdatCodes);
                 else
                     FD.xdat(lastxdn + k).id = xdatmap(FDin(i).xdat(k).id);                     
                     FD.xdat(lastxdn + k).code = FD.xdatCodes{xdatmap(FDin(i).xdat(k).id)};                     
                 end
           end
           
          if isfield(FD,'sac')
             for k = 1:length(FDin(i).sac)

%                  FD.sac(lastsacn + k) = FDin(i).sac(k);
                 FD.sac(lastsacn + k).startT = double(FDin(i).sac(k).startT) + lastendt;
%                  FD.sac(lastsacn + k).endT = FDin(i).sac(k).endT + lastendt;
                 FD.sac(lastsacn + k).xdhist = xdatmap(FDin(i).sac(k).xdhist);
                 FD.sac(lastfxn + k).xdindex = FDin(i).sac(k).xdindex + lastxdn;
             end    
          end
          
          
%         end

        lastendt = double(FD.fix(lastfxn + k).endT) + tpad ;
        lastendts(i) = double(lastendt);
        lastfxn = fixns(i)+lastfxn;
        lastfxns(i) = lastfxn;
        lastsacn = sacns(i)+lastsacn;
        lastxdn = xdns(i) + lastxdn;
        if isnumeric(FDin(i).xdatCodes)
            xdmaps{i}(FDin(i).xdatCodes) = FD.xdatCodes(xdatmap);
        else
            xdmaps{i} = FD.xdatCodes(xdatmap);
        end
end

        