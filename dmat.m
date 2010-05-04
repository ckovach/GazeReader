function [D,Drank] = dmat(X,r)

%function [D,Drank] = dmat(X)
% Computest the Euclidean distance between all rows of X and returns both the
% distance matrix, and the distance rank matrix;
%
%function [D,Drank] = dmat(X,r)
%   Uses Lr norm with parameter r. r=1 is city block distance, r=2 is Euclidean
%

%C. Kovach 2008

if nargin < 2
    r = 2;
end

nrows = size(X,1);
X1 = repmat(permute(X,[1 3 2]),1,nrows);
Xt = repmat(permute(X,[3 1 2]),nrows,1);

D = sum((X1 - Xt).^r,3).^(1/r);

[srt,Drank] = sort(D);


