function ptch = pcolor_radial(TH,R,C,ax)

%Makes radial surface plot

if nargin < 4
    ax = gca;
end

[th,r] = meshgrid(TH,R);

vtx = kron(ones(size(th,2)-1,1),[1:size(th,1)-1]') + size(th,1)*(kron([1:size(th,2)-1]',ones(size(th,1)-1,1))-1);

vertices = [vtx,vtx+1,vtx+size(th,1)+1,vtx+size(th,1)]';

X = cos(th(:)).*r(:);
Y = sin(th(:)).*r(:);

axes(ax)
ptch = patch(X(vertices),Y(vertices),C(vertices));

shading flat