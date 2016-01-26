% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


function Xtrans = pdf_point_transform(PDF,X)


if length(PDF) > 2 || min(size(PDF)) > 1
    
    prx = cumsum(sum(PDF));
    pry = cumsum(sum(PDF,2))';
else
    prx = cumsum(PDF);
    pry = 1;
end
nx = length(prx);
ny = length(pry);
npt = size(X,1);

dx = repmat(X(:,1),1,nx) - repmat(prx,npt,1);
dy = repmat(X(:,2),1,ny) - repmat(pry,npt,1);


[mnx,mnxi] = min(abs(dx),[],2);    
[mny,mnyi] = min(abs(dy),[],2);

Xtrans = [mnxi,mnyi];

