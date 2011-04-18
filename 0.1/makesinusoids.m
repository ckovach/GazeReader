
function S = makesinusoids(X,order)

% Creates a set of 2*order sinusoids of increasing frequency and period 1.
% X may be a column vector or N x 2 matrix. Frequencies are bounded at
% pi*order.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


[frqfun,fr] = makefourier(order);

S = frqfun(X,1);
