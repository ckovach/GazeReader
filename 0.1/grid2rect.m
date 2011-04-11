function rect = grid2rect(grid)

% rect = grid2rect(grid)
%   Converts from 2 element cell array specification of a bin grid to an
%   (M*N)x4 rect matrix. grid{1} is the position of the grid in rect format 
%   and grid{2} is [N M], the number of divisions along the X and Y axes respectively.
% 
% See also RECT2GRID

if ~iscell(grid)
    rect = grid;
    return
end


lx = linspace(grid{1}(1), grid{1}(2), grid{2}(1)+1);
lxoffset = (grid{1}(2) - lx(end))./2;
lx = lx+lxoffset;

ly = linspace(grid{1}(3), grid{1}(4), grid{2}(2)+1);
lyoffset = (grid{1}(4) - ly(end))./2;
ly = ly+lyoffset;

[X,Y] = meshgrid(lx(1:end-1),ly(1:end-1));

[dx,dy] = meshgrid(diff(lx),diff(ly));

rect = cat(2 , X(:) , X(:) + dx(:) ,  Y(:), Y(:) + dy(:) );

