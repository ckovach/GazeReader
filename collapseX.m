
function [ I, nrows, uindx, varargout] = collapseX(varargin)


% [Xunq, nrows, I] = collapseX(X)
% 
% Returns unique rows of X, frequency of occurrence, and index, 
% where Xunq = X(I).
%

% C. Kovach 2011

X = [varargin{:}];

    
[Xunq,I,uindx] = unique(X,'rows'); 
 
nrows = zeros(size(Xunq,1),1);

dsrtu = diff([sort(uindx)]);
nrows = diff([0;find(dsrtu);length(dsrtu)]);

if nargout > 2
    for i = 1:length(varargin)
        varargout{i} = varargin{i}(I,:);
    end
end
    
    
    
    