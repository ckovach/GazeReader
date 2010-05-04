function ED = concatED(lastendts,xdmaps,varargin)

%EDcat = concatED(lastendts,xdmaps, ED1,ED2,ED3,...)
%Concatenates multiple eventsationData structures into a single one.


if length(varargin)>1
    EDin = cat(2,varargin{:});
else
    EDin = varargin{1};
end

lastendt = 0;
lastcode = 0;
lastevn = 0;
lastxdn = 0;

ED = EDin(1);
ED.events = cat(2,EDin.events);
ED.xdat = cat(2,EDin.xdat);

for i = 1:length(EDin)

        if i == 1
            ED= EDin(i);
        else

            
             for k = 1:length(EDin(i).events)   

%                  ED.events(end + k) = EDin(i).events(k);
                 ED.events(lastevn + k).time = EDin(i).events(k).time + lastendt;
                 ED.events(lastevn + k).code = EDin(i).events(k).code + lastcode;
                 ED.events(lastevn + k).xdatcode= xdmaps{i}(ismember(xdmaps{i},EDin(i).events(k).xdatcode));
             end
          
           for k = 1:length(EDin(i).xdat)   
%                  ED.xdat(end + k) = EDin(i).xdat(k);
                 ED.xdat(lastxdn + k).startT = EDin(i).xdat(k).startT + lastendt;
                 ED.xdat(lastxdn + k).id = xdmaps{i}(ismember(xdmaps{i}, EDin(i).xdat(k).id) );
           end
           
          
          
        end

        lastendt = lastendts(i);
end

ED.codeincr = ED.events(end).code;



        