function tuntil = time_until(delta_train)

% Takes as input a train of delta functions and returns a vector with the time
% lag until the following delta in sampling units

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


delt = find(delta_train);

xt = zeros(size(delta_train));
xt([1;delt(1:end-1)]) = delt - [0;delt(1:end-1)];

tuntil = cumsum(xt)-repmat((1:size(delta_train,1))',1,size(xt,2));

