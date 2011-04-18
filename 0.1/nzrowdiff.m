
function q = nzrowdiff(X)

% Returns a matrix the size of X in which each nonzero element is replaced
% by the number of rows between it and the nearest preceding non-zero
% element within the column. The first non-zero element in each column is
% replaced with zero. Eg.
% 
% >> X = [1 0 3 0; 0 0 6 -1; 0 2 0 1; 7 0 5 0]
% 
% X =
% 
%      1     0     3     0
%      0     0     6    -1
%      0     2     0     1
%      7     0     5     0
% 
% >> Y = nzrowdiff(X)
% 
% Y =
% 
%      0     0     0     0
%      0     0     1     0
%      0     0     0     1
%      3     0     2     0

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


fm = sparse(X~=0);
q = zeros(size(fm));
[ir,jc] = sparseij(fm);
di = [0;diff(ir)];
di(jc(1:end-1)+1) = 0;
q(fm) = di;