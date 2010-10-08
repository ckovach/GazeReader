
function R = makeregressor(X,varargin)

%  makeregressor(X,'noptions',nopt)
%
% 
%   MAKEREGRESSOR takes input X and puts it into the form of a regressor
%   structure.
%   
%   Nopt is the the number-of-options vector giving the number of options
%   (bins) for each event. X must have sum(Nopt) rows if Nopt is a vector
%   and have an integer multiple of Nopt if Nopt is a scalar.
%
%   For w = cumsum([0,Nopt]), the row in X associated with the jth option on
%   the ith trial is given by w(i) + j if Nopt is a vector, and
%   Nopt*(i-1)+j if Nopt is scalar.
%
%   makeregressor([])
%     
%   Initializes an empty regressor
%
% See also INTERACTION, SPLIT, POOL

%
% C. Kovach 2007
% 

sparseform = false;
noptions = [];
fixed =[];
interaction_order = 1;
i=1;
codeincr = 0;
normconst = 1;
label = '';
mkfun = @(X)X;
noptions = [];
Npar = 0;
postmult = ones(size(X,1),1);
while i <= length(varargin)
    
    switch lower(varargin{i})
        
        case 'sparseform' %Value is represented as a sparse  block diagonal matrix 
            sparseform = 1;
        case 'fullform' %Value is represented as a full matrix with blocks catenated along the 2nd dimension
            sparseform = 0;
        case 'normconst' %Value is represented as a sparse  matrix if true
            normconst = varargin{i+1};
            i=i+1;
        case 'codeincr' 
            codeincr = varargin{i+1};
            i=i+1;
        case 'noptions' %Vector of number of outcomes possible on each trial (or scalar constant)
            noptions = varargin{i+1};
            i=i+1;
        case 'label' %Vector of number of outcomes possible on each trial (or scalar constant)
            label = varargin{i+1};
            i=i+1;
        case 'fixed' %Vector of parameter estimates which are not to be fitted. All values = zero are fitted (to fix at zero let value = eps)
            fixed = varargin{i+1};
            i=i+1;
        case 'interactionorder' %If this regressor represents an interaction, specify its order. 
            interaction_order = varargin{i+1};
            i=i+1;
        case 'function' %Applies specified function the input
            mkfun = varargin{i+1};
            i=i+1;
        case 'postmultiply' %after the regressor matrix is computed multiply each column by this vector
            postmult  = varargin{i+1};
            i = i+1; 
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end

if nargin < 2
    label = [];
end

R = struct('value',[],'info',[],'normconst',normconst,'label',label,'noptions',...
            noptions,'Npar',Npar,'code',codeincr+1,'codevec',[],'factmat',[],'levmat',[],...
             'fixed',[],'function',[],'deriv',[]);
R.info = struct('form',[],'intxnord',[],'COMMAND',[],'hashcode',[],'parent',[],'label',...
                [],'pooled_labels',{{}},'contrasts',{{}},'parentindex',[]);
R.Npar = 0;
if isempty(X) %Initializes an empty regressor if X is empty
    return
end

% X = X.*normconst;
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

MKX = mkfun(X);

if length(postmult) > 1
    MKX = MKX.*repmat(postmult,1,size(MKX,2));
else
    MKX = MKX.*postmult;
end

npar = size(MKX,2);


if sparseform && ~issparse(X)
    
    R.value = sparseblock(MKX + eps ,noptions,'transpose');
    R.info.form = 'sparse';
else
    
    R.value = MKX + eps;
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

R.function = mkfun;


reseed;

R.info.hashcode = uint32(rand*intmax('uint32')); %For now a random number will be the hashcode

R.code = codeincr + 1;

R.factmat = ones(1,npar)*R.code;
R.levmat  = 1:npar;

R.codevec = ones(1,npar)*R.code;

if isempty(fixed) 
    fixed = zeros(1,R.Npar);
elseif length(fixed)~=R.Npar
    error('Fixed vector must be same length as number of parameters or empty');
end

R.fixed = fixed;

% nzv = typecast(nonzeros(R.value)-mean(nonzeros(R.value)),'int32');
% rand('seed',234952);
% sumvec = sign(rand(nnz(R.value),1)-.5);
% R.info.hashcode = sum(  nzv(sumvec>0,:), 'native' )./2 + sum( nzv(sumvec<0,:), 'native' )./2;


