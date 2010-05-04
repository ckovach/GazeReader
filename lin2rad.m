function [R,Th] = lin2rad(X,Y)

%Converts linear XY coordinates to polar coordinates with Th = 0 being
%horizontal up.

if nargin == 2
    X = [X(:),Y(:)];
end

R = sqrt(sum(X.^2,2));

Th = atan2(X(:,1),X(:,2));