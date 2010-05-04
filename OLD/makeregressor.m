
function R = makeregressor(X,label,varargin)

% function makeregressor(X,label,normconst)
% 
%   MAKEREGRESSOR takes input X and puts it into the form of a regressor
%   structure.
% 
% function makeregressor([])
%     
%   Initializes an empty regressor
%
% See also INTERACTION, SPLIT, POOL

%
% C. Kovach 2007
% 

normconst = 1;

sparseform = 1;
noptions = [];
interaction_order = 1;
i=1;
while i <= length(varargin)
    
    switch lower(varargin{i})
        
        case 'sparseform' %Value is represented as a sparse  block diagonal matrix 
            sparseform = 1;
        case 'fullform' %Value is represented as a full matrix with blocks catenated along the 2nd dimension
            sparseform = 0;
        case 'normconst' %Value is represented as a sparse  matrix if true
            normconst = varargin{i+1};
            i=i+1;
        case 'noptions' %Vector of number of outcomes possible on each trial (or scalar constant)
            noptions = varargin{i+1};
            i=i+1;
        case 'interactionorder' %If this regressor represents an interaction, specify its order. 
            interaction_order = varargin{i+1};
            i=i+1;
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end

if nargin < 2
    label = [];
end

R = struct('value',[],'info',[],'normconst',[],'label',[],'noptions',[],'Npar',[]);
R.info = struct('form',[],'intxnord',[],'COMMAND',[],'hashcode',[],'parent',[],'label',[],'pooled_labels',{{}},'contrasts',{{}});

if isempty(X) %Initializes an empty regressor if X is empty
    return
end
if size(X,3) > 1
    noptions = size(X,2);   
    X = reshape(X,[size(X,1), size(X,2)*size(X,3)])';
end


if length(noptions) == 1
    
        ntrials = size(X,1)./noptions;
        
    if  ntrials - round(ntrials)~= 0
        error('noptions is not a divisor of size(X,1)')
    end
    
    noptions = noptions*ones(ntrials,1);
end

ntrials = length(noptions);
npar = size(X,2);


if sparseform && ~issparse(X)
    
    R.value = sparseblock(X + eps ,noptions,'transpose');
    R.info.form = 'sparse';
else
    
    R.value = X + eps;
    R.info.form = '2d_block';

end

R.normconst = normconst;

R.label = label;
R.info.label = label;

R.noptions = noptions;
R.Npar = npar;
R.info.intxnord = interaction_order;
try
    fid = fopen([mfilename,'.m'],'r');
    R.info.COMMAND =  fread(fid);
    fclose(fid);
catch
    warning('Failed to record the script used to generate regressor %s poly',R.label);
end


reseed;

R.info.hashcode = uint32(rand*intmax('uint32')); %For now a random number will be the hashcode

% nzv = typecast(nonzeros(R.value)-mean(nonzeros(R.value)),'int32');
% rand('seed',234952);
% sumvec = sign(rand(nnz(R.value),1)-.5);
% R.info.hashcode = sum(  nzv(sumvec>0,:), 'native' )./2 + sum( nzv(sumvec<0,:), 'native' )./2;

