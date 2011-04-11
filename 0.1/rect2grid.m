function grid = rect2grid(rect)

% grid = rect2grid(rect)
%   From an (M*N)x4 matrix specifying a grid, rect, generates a cell array, grid, where grid{1}
%   is the position of the grid in rect format and grid{2} is [N M], the
%   number of divisions along the X and Y axes respectively.
% 
% See also GRID2RECT

if iscell(rect) || isempty(rect)
    grid = rect;
    return
end

minmax = cat(1,min(rect,[],1),max(rect,[],1));

grid{1} = minmax([1 4 5 8]);

nbins = [length(unique(rect(:,1))),length(unique(rect(:,3)))];
grid{2} = nbins;

