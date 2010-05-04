
function R = buildpolyreg(x,polyord,varargin)

% R = buildpolyreg(x,polyord)
% 
%   Generates a regressor block of polynomial terms in x, where x is a N x 1 vector
%
%   Each block is is [x,x^2,x^3,...,x^polyord]
% 
%   


normconst = 1; 
label = '';

i = 1;
minord = 1;
polyord_is_vector = false;

trig = false;

while i <= length(varargin)
    
    switch lower(varargin{i})
        
        case 'normconst' %multiply by this value before building polynomials.
            normconst = varargin{i+1};
            i = i+1;
        case 'label' %prepend this string to the regressor label
            label  = varargin{i+1};
            i = i+1; 
        case 'trig' %Polynomials are trigonometric, that is, sinusoids.
            trig = true;
        case 'dc' %include the "DC" zeroth order term (default is to start at 1st order)
            minord  = 0;
        case 'polyvec' %treat polyord as a 1xp vector, explicitly representing terms in the polynomial
                    %(for building polynomials that skip some terms)
               polyord_is_vector = true;
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end

if polyord_is_vector
    poly = polyord;
else
    
    poly = minord:polyord;
end

if isempty(label)
    label = 'Poly';
else
    label = strcat(label,'_poly');
end

R = makeregressor(repmat(x,1,length(poly)).^repmat(poly,size(X,1),1),label,'normconst',normconst);

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
    fid = fopen(mfilename,'r');
    R.info.COMMAND =  fread(fid);
    fclose(fid);
catch
    warning('Failed to record the script used to generate regressor %s poly',R.label);
end


