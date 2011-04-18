function tsince= time_since(delta_train)

% Takes as input a train of delta functions and returns a vector with the time
% since the preceding delta in sampling units


delt = find(delta_train);

xt = zeros(size(delta_train,1),1);
xt([delt(1:end)]) = delt - [0;delt(1:end-1)];

tsince= (1:size(delta_train,1))'-cumsum(xt);
tsince(1:delt(1)-1) = 0;


