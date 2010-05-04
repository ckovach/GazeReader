
function R = buildbinreg(x,b,varargin)

% R = buildbinreg(x,b)
% 
%   Generates a regressor block of bins over x, where x is a N x 1 vector
%
%   Column i has the value x > b(i-1) & x < b(i), where b(0) = -Inf and
%   b(end+1) = Inf.
% 

normconst = 1; 
label = '';

i = 1;
% minord = 1;
% polyord_is_vector = false;
% 
% trig = false;
noptions = [];
codeincr = 0;
subtractdc = false;
postmult = ones(size(x,1),1);
center  = false;
binN = 0; %When non-zero, then bin boundaries are adjusted so that each bin contains binN objects.
baseline = 0;
if length(b) == 1
    binN = b;
end

while i <= length(varargin)
    
    switch lower(varargin{i})
        
        case 'normconst' %multiply by this value before building polynomials.
            normconst = varargin{i+1};
            i = i+1;
        case 'label' %prepend this string to the regressor label
            label  = varargin{i+1};
            i = i+1; 
        case 'noptions' %number of options vector
            noptions  = varargin{i+1};
            i = i+1; 
        case 'postmultiply' %after the regressor matrix is computed multiply each column by this vector
            postmult  = varargin{i+1};
            i = i+1; 
        case 'codeincr' %prepend this string to the regressor label
           codeincr = varargin{i+1};
            i = i+1; 
        case 'binn' %select binning so that each bin contains this many objects
           binN = varargin{i+1};
            i = i+1; 
        case 'center' %Subtract mean value for each bin
            center  = true;
        case 'subtractdc' %Subtract mean value for each bin
            subtractdc  = true;
        case 'baseline' %Baseline
           baseline =  varargin{i+1};
            i = i+1; 
       
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end

% persistent poly frqfun

if isempty(label)
    label = 'Bins';
end


if binN ~=0
    srtx = sort(x);

    b = srtx(binN:binN:end-binN)';
    db = diff(b);
    db(1) = 1;
    b(db==0) = [];
end

if baseline > 0
    zmat = eye(length(b)+1);
    zmat(:,baseline) = [];
    polyfun = @(x) (repmat(x,1,length(b)+1) < repmat([b,Inf],size(x,1),1) & repmat(x,1,length(b)+1) > repmat([-Inf, b],size(x,1),1))...
                    *zmat;
    ncol = length(b);
elseif subtractdc
    db = diff([0,b,1]);
    polyfun = @(x) (repmat(x,1,length(b)+1) < repmat([b,Inf],size(x,1),1) & repmat(x,1,length(b)+1) > repmat([-Inf, b],size(x,1),1))...
                - repmat(db,size(x,1),1);
            
    ncol = length(b)+1;
else
    polyfun = @(x) repmat(x,1,length(b)+1) < repmat([b,Inf],size(x,1),1) & repmat(x,1,length(b)+1) > repmat([-Inf, b],size(x,1),1);
    ncol = length(b)+1;
end

if center
    regfun = @(x) polyfun(x) - ones(size(x,1),1)*mean(polyfun(x),1);
else
    regfun = polyfun;
end



R = makeregressor(regfun(x).*repmat(postmult,1,ncol),'label',label,'normconst',normconst,'noptions',noptions,'codeincr',codeincr);

R.function = @(X) makeme(X,regfun,b);

try
    fid = fopen([mfilename,'.m'],'r');
    R.info.COMMAND =  fread(fid);
    fclose(fid);
catch
    warning('Failed to record the script used to generate regressor %s poly',R.label);
end

%--------------------------
function R = makeme(X,regfun,b)

% if size(X,2) < length(b)+1   
%     error('Input matrix must have %i columns',length(b)+1)
% end

R = regfun(X);

