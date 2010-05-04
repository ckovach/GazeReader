

function [FixBin,EQ,BinPos,gridedges] = MakeGrid( FixPos, gridsize, gridRange)

%Creates a grid of specified size and identifies bin for each observation


if nargin < 2 || isempty(gridsize)
    gridsize = [16 16];
elseif length(gridsize) == 1
    gridsize = gridsize*[1 1];
end    
    

if nargin < 3
    % range covered by the fine grid [xmin xmax ymin ymax]
    % gridRange = [-.25 1.25  -.25 1.25];
    gridRange = [-0 1 0 1];
end
dgr = diff(gridRange);
dgr = dgr([1,3]);

gridedges = {gridRange(1)+dgr(1)/gridsize(1):dgr(1)/gridsize(1):gridRange(2),gridRange(3)+dgr(2)/gridsize(2):dgr(2)/gridsize(2):gridRange(4)} ;
%Bin Positions
[BinX , BinY] = meshgrid( gridedges{1},gridedges{2});
BinX = BinX - .5*dgr(1)/gridsize(1);
BinY = BinY - .5*dgr(2)/gridsize(2);
BinPos = [BinX(:),BinY(:)];

gridedges = {gridRange(1):dgr(1)/gridsize(1):gridRange(2),gridRange(3):dgr(2)/gridsize(2):gridRange(4)} ;

%Bin associated with each fixation
FixBin = floor( ( FixPos(:,1)- gridRange(1))./dgr(1) * gridsize(1) )*gridsize(2) + ceil(( FixPos(:,2) - gridRange(3))./dgr(2) * gridsize(2) );

%Adding bins corresponding to external quadrants
% Left Upper Quadrant outside of face
EQ( (FixPos(:,1) < gridRange(1) & FixPos(:,2) < .5 ) | (FixPos(:,1) < .5 & FixPos(:,2) < gridRange(3) ) ) = length(BinPos )+1;
% Left Lower Quadrant outside of face
EQ( (FixPos(:,1) < gridRange(1) & FixPos(:,2) >= .5 ) | (FixPos(:,1) < .5 & FixPos(:,2) >= gridRange(4) ) ) = length(BinPos )+2;
% Rigth Lower Quadrant outside of face
EQ( (FixPos(:,1) > gridRange(2) & FixPos(:,2) >= .5 ) | (FixPos(:,1) > .5 & FixPos(:,2) >= gridRange(4) ) ) = length(BinPos )+3;
% Right Upper Quadrant outside of face
EQ( (FixPos(:,1) > gridRange(2) & FixPos(:,2) < .5 ) | (FixPos(:,1) > .5 & FixPos(:,2) < gridRange(3) ) ) = length(BinPos )+4;

% BinPos(end+1:end+4,:) = 0; %Ignoring position for exterior bins

