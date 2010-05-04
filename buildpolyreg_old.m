
function R = buildpolyreg(x,polyord,varargin)

% R = buildpolyreg(x,polyord)
% 
%   Generates a regressor block of polynomial terms in x, where x is a N x 1 vector
%
%   Each block is is [x,x^2,x^3,...,x^polyord]
% 
% R = buildpolyreg(x,polyord,'trig')
%   
%   Returns trigonometric polynomials (sinusoids)
% 
%   R = [sin(x), cos(x), sin(2*x), cos(2*x),...,sin(nx),cos(nx)]


normconst = 1; 
label = '';

i = 1;
minord = 1;
polyord_is_vector = false;

trig = false;
noptions = [];
codeincr = 0;
% subtractdc = false;
subtractdc = true;
postmult = ones(size(x,1),1);
center  = false;

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
        case 'subtractdc' %subtracts a constant term equivalent to decorrelating with dc
            subtractdc = varargin{i+1};
            i = i+1; 
        case 'codeincr' %prepend this string to the regressor label
           codeincr = varargin{i+1};
            i = i+1; 
        case 'trig' %Polynomials are trigonometric, that is, sinusoids.
            trig = true;
        case 'dc' %include the "DC" zeroth order term (default is to start at 1st order)
            minord  = 0;
            subtractdc = false;
        case 'center' %Subtract mean value for each polynomial term except DC (if it's included)
            center  = true;

        case 'polyvec' %treat polyord as a 1xp vector, explicitly representing terms in the polynomial
                    %(for building polynomials that skip some terms)
               polyord_is_vector = true;
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end

% persistent poly frqfun

if polyord_is_vector
    poly = polyord;
else
    
    poly = minord:polyord;
end

if isempty(label)
    label = 'Poly';
% else
%     label = strcat(label,'_poly');
end

if subtractdc 
    subdc = 1./(poly + 1);
else
    subdc = zeros(1,length(poly)); 
end


if ~trig    
    polyfun = @(x) (repmat(x,1,length(poly)).^repmat(poly,size(x,1),1) - repmat(subdc,length(x),1)); %currently this only creates 1d polynomials
    ncol = length(poly);
else
     mkfarg = {};
    if ~subtractdc 
        mkfarg = {'keepdc'};
    end
    [makefourier_output,fr] = makefourier(polyord*ones(1, size(x,2)),mkfarg{:});     %trig polynomials may be in an arbitrary number of dimension based on length of polyord
    polyfun = @(x)makefourier_output(x,1)';    
%     ncol = 2*length(poly)-sum(poly==0);
    ncol = length(polyfun(x(1,:)));
end    

if center
    regfun = @(x) polyfun(x) - ones(size(x,1),1)*mean(polyfun(x),1)*diag((minord:polyord)>0);
else
    regfun = polyfun;
end

R = makeregressor(regfun(x).*repmat(postmult,1,ncol),'label',label,'normconst',normconst,'noptions',noptions,'codeincr',codeincr);

R.function = @(X) makeme(X,regfun,polyord);

% R.value = repmat(x,1,length(poly)).^repmat(poly,size(X,1),1);
% R.normconst = normconst;
% 
% if isempty(label)
%     R.label = 'Poly';
% else
%     R.label = strcat(label,'_poly');
% end
% 
% R.info.polyterms = poly;

try
    fid = fopen([mfilename,'.m'],'r');
    R.info.COMMAND =  fread(fid);
    fclose(fid);
catch
    warning('Failed to record the script used to generate regressor %s poly',R.label);
end

%--------------------------
function R = makeme(X,regfun,polyord)

if size(X,2) < length(polyord)   
    error('Input matrix must have %i columns',length(polyord))
end

R = regfun(X);

