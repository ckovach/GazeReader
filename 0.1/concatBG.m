function [BG,codemap] = concatBG(varargin)

%BGcat = concatBG(BG1,BG2,BG3,...)
%Concatenates multiple binData structures into a single one.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


if length(varargin)>1
    BGin = cat(2,varargin{:});
else
    BGin = varargin{1};
end



lastcode = 0;
for i = 1:length(BGin)

        if i == 1
            BG= BGin(i);
        else

             for k = 1:length(BGin(i).groups)   

                 BG.groups(end + k) = BGin(i).groups(k);
                 BG.groups(end).code = BGin(i).groups(k).code + lastcode;
             end
           
          
        end
        if ~isempty(BGin(i).groups)
            oldcodes = [BGin(i).groups.code];   
        end 
           
        newcodes = oldcodes+lastcode;
        codemap{i}(oldcodes) = newcodes;

        lastcode = BG.groups(end).code;
        
end

BG.codeincr = BG.groups(end).code;
        