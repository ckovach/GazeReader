
function [dbingr1,dbingr2] = bindifference(bingr1,bingr2)

% function [ubigr1,uibgr2] = bindifference(bingr1,bingr2)
% 
% Returns 1 for bins in  bingr 1 whose position variable falls within no members
% of bingr 2, and 1 for bins in bingr2 which contain no bins in group 1. 
%
%
% See also MAKEBINDATA and BINUNION

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

for j = 1:length(bingr2)
    dbingr2{j} =1;
end

for i = 1: length(bingr1)
    dbingr1{i} = 1;
    for j = 1:length(bingr2)
        ismem = bingr2(j).isinside(bingr2(j),bingr1(i).centers);

        dbingr1{i} = (sum(ismem,2) == 0) & dbingr1{i};

        dbingr2{j} = sum(ismem,1) == 0 & dbingr2{j};
    end
end
