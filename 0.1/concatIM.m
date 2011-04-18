function [IM,codemap] = concatIM(varargin)

%IMcat = concatIM(IM1,IM2,IM3,...)
%Concatenates multiple imageData structures into a single one.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


if length(varargin)>1
    IMin = cat(2,varargin{:});
else
    IMin = varargin{1};
end



lastcode = 0;
for i = 1:length(IMin)

        if i == 1
            IM= IMin(i);
        else

             for k = 1:length(IMin(i).images)   

                 IM.images(end + k) = IMin(i).images(k);
                 IM.images(end).code = IMin(i).images(k).code + lastcode;
             end
           
          
        end
        if ~isempty(IMin(i).images)
            oldcodes = [IMin(i).images.code];
        end
        newcodes = oldcodes+lastcode;
        codemap{i}(oldcodes) = newcodes;

        lastcode = IM.images(end).code;
        
end

IM.codeincr = IM.images(end).code;
        