function [Rout,interleaver] = interleave(varargin)

% [Rout,interleaver] = interleave(R1,R2,...)
%
% Interleaves the rows of the design matrices in regressor structures
% R1,R2,.... The output is a regressor for which outcomes are are concatenated
% so that Rout.noptions(k) = R1.noptions(k) + R2.noptions(k) +  ... 
%
% See also MAKEREGRESSOR, POOL, SPLIT




Rs = [varargin{:}];



noptmtx = [Rs.noptions];
noptions = sum(noptmtx,2);

spe = speye(sum(noptions)); 

X = 0;

index = (1:sum(noptions))';

ci = zeros(sum(noptions),1);
ci(cumsum(noptions(1:end-1))+1) = noptions(1:end-1);

expander = cumsum([1;ci>0]);
expander(end) = [];



rowindex = index-cumsum(ci);


cnopts = cumsum([zeros(size(noptmtx,1),1),noptmtx],2);


npar = Rs(1).Npar;

X = zeros(sum(noptions),npar);


interleaver = @(varargin) interleavfun(cnopts,expander,rowindex,varargin{:});

X = interleaver(Rs.value);
lbl = sprintf('%s > ',Rs.label);
lbl(end) = [];

Rout = makeregressor(X,'noptions',noptions,'label',lbl);

% for i = 1:length(Rs)
%     
%     X((cnopts(expander,i) < rowindex) & (rowindex < cnopts(expander,i+1)+1),:) = Rs(i).value;
%     
% end    
    

%%%%

function X = interleavfun(cnopts,expander,rowindex,varargin)
    
X = zeros(size(expander,1),size(varargin{1},2));

for i = 1:length(varargin)
    
    X((cnopts(expander,i) < rowindex) & (rowindex < cnopts(expander,i+1)+1),:) = varargin{i};
    
end    
    