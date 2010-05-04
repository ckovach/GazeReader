
function [ubingr1,ubingr2] = binunion(bingr1,bingr2)

% function [ubigr1,uibgr2] = binunion(bingr1,bingr2)
% 
% Returns 1 for bin in  bingr 1 whose center falls within members
% of bingr 2, and 1 for bins in bingr2 which contain a bin in group 1. 
%
% See also MAKEBINDATA and BINDIFFERENCE

for i = 1: length(bingr1)
    dbingr1{i} = 0;
    for j = 1:length(bingr2)
        ismem = bingr2(j).isinside(bingr2(j),bingr1(i).centers);

        ubingr1{i} = sum(ismem,2) > 0 | dbingr1{i};

        ubingr2{j} = sum(ismem,1) > 0;
    end
end


