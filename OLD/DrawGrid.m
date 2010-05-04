function DrawGrid(gridedges,varargin)


xlim = gridedges{1}([1 end]);
ylim = gridedges{2}([1 end]);
hold on,
for i = 1:length(gridedges{1}) 
    plot(gridedges{1}(i)*[1 1],ylim,varargin{:});
end

for i = 1:length(gridedges{2}) 
    plot(ylim,gridedges{2}(i)*[1 1],varargin{:});
end
   