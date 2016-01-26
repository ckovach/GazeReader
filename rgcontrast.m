
function C = rgcontrast(R,code)

%Returns an n x k contrast matrix that gives parameter estimates
%for the regressor group specified by code.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


allcodes = [R.codevec];

C = diag(ismember(allcodes,code));
C = C(:,sum(C)>0);

