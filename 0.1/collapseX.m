
function [ Iout, nrows, uindxout,bout, varargout] = collapseX(X,b,varargin)


% [I, nrows, index,noptsout] = collapseX(X,noptsin,X2,X3,...)
% 
% Returns and index into unique row blocks of [X,X1,X2,...], frequency of occurrence (count),  
% index of rows into Xunq, where Xunq = X(I,:) and X = Xunq(index,:), and noptions vector. 
%
% Input noptsin contains the number of rows in each succesive block.
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2011



moreX = {};
% discardtrials = false;
i = 1;
discard = false;
while i <= length(varargin)
    
    if ~ischar(varargin{i})
        moreX = [moreX,varargin{i}];
        
    elseif ischar(varargin{i})
        switch lower(varargin{i})
            case 'discardtrials'
                discard = varargin{i+1};
                i = i+1;          
            otherwise
                error(sprintf('%s is not a valid keyword',varargin{i}));
        end

    end
    
       i = i+1;
    
end



if nargin < 2
    b = 1;
end

X = cat(2,X,moreX{:});




if length(size(b)) > 2 || ~any(size(b) == 1)
    error('Second argument must be a vector.');
end
    
b = b(:);
nblocks = length(b);

if nblocks == 1
    nblocks = size(X,1)./b;
    if mod(nblocks,1) ~= 0
        error('b doesn''t evenly divide the number of Columns of X.')
    end
    b = b*ones(nblocks,1);
end

if sum(b)~=size(X,1)
    error('Total bin count doesn''t match the number of rows in X.')
end


if max(b) > 1   %For blocks with multiple rows in X, we must identify unique rows
      
    
    
   
    %Index arrays so we can properly reorder things
    trindex = zeros(size(X,1),1);
    trindex([1;cumsum(b(1:end-1))+1]) = 1;
    trindex = cumsum(trindex);

    %Discard trials for which noptions is zero
    if any(b==0)
        discard = discard | b==0;
    end
    
%     if any(discard)
%         X(discard(trindex) ,:) = [];
%         trindex(discard(trindex),:) = [];
%         b(discard) = [];
%     end

    csb = cumsum(b(1:end-1));
    blindex = ones(size(X,1),1);
    blindex(csb+1) = blindex(csb+1)-b(1:end-1);
    blindex = cumsum(blindex);
    
    csb = [0;csb];
    
    xindex0 = csb(trindex);

    
    %%% Using some sparse matrix acrobatics we can efficiently reorder the
    %%% blocks into rows.
    
    xtest = X';
    
    xspb = sparseblock(xtest(:),b*size(X,2),'transpose');
    [ir,jc] = sparseij(xspb);
    jind = zeros(size(ir));
    jind(jc(2:end)+1) = diff(jc);
    csj = cumsum(jind);
    irnew = ir-csj(1:end-1);
    nrows = max(irnew)+1;
    Xrow = sparseblockmex(xtest(:),irnew,jc, nrows)'; % Blocks reordered into zero-padded rows.

    
else
    
    Xrow = X;
end

if any(discard)
    Xrow(discard,:) = 0;
    Xrow = [discard,Xrow];
end
    
[~,I,uindx] = unique(Xrow,'rows'); 

if any(discard)
    I(discard(I)) = [];
    uindx(discard) = [];
end


dsrtu = diff([sort(uindx)]);
nrows = diff([0;find(dsrtu);length(dsrtu)+1]);

bout = b(I);

if max(b)>1
   %%% Adjust indices to the original structure of X.
   csbi = csb(I);   

    csbout = cumsum(bout(1:end-1));
    blindexi = ones(sum(bout),1);
    blindexi(csbout+1) = blindexi(csbout+1)-bout(1:end-1);
    blindexi = cumsum(blindexi);
    
   expindx = sparseblock(ones(sum(bout),1),bout,'transpose');
   
   
   Iout = expindx*csbi + blindexi;
   
   uindxout = uindx(trindex);
   
else
    
    Iout = I;
    uindxout = uindx;
    
end


if nargout > 4
    outargs = [{X},moreX{:}];
    for i = 1:length( outargs)
        varargout{i} = outargs{i}(Iout,:);
    end
end
    
    
    
    