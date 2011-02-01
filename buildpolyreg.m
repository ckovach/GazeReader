
function R = buildpolyreg(x,polyord,varargin)

% R = buildpolyreg(x,polyord)
% 
%   x is an Nxk input matrix and polyord is a 1xk matrix of polynomial
%   orders.
%
%   For k == 1
%   R.value = [x,x^2,x^3,...,x^polyord]
% 
%   For k > 1
%   R.value = [..., x(:,1)^a1*x(:,2)^a2*x(:,3)^a3..., ...]
%   for all values of a1,a2,a3,...ak such that ai < polyord(i), and 
%    a1+a2+...+ak <= max(polyord).
% 
% 
% R = buildpolyreg(x,polyord,'trig')
%   
%   Return series of 1D or 2D trigonometric polynomials (that is, sinusoids)
% 
%   R.value = [sin(x), cos(x), sin(2*x), cos(2*x),...,sin(nx),cos(nx)]
%
% R = buildpolyreg(x,polyord,'option',value)
%
% Create the regressor with any of the following options:
%      'premultiply' -  multiply the columns of the input matrix by this value
%                         BEFORE constructing the polynomial.
%       'postmultiply' -  multipy the columns of the output matrix AFTER 
%                           constructing the polynomial.
%       'subtractdc'   - subtracts a constant term equivalent to
%                        decorrelating with dc over the unit invterval (default = true). 
%       'trig'         -  No argument. Polynomials are trigonometric, that is, sinusoids.
%       'dc'           - No argument. Include the unit constant (zeroth order) term 
%                        (default is to start at 1st order)
%       'center'      - Subtract mean value for each polynomial term except DC (if it's included)
%       'clip'        - for values of the input outside this range, the row of the output matrix is set to zero
%       'polyvec'     -  treat polyord as a 1xp vector, explicitly representing 
%                           terms in the polynomial (for building
%                           polynomials that skip some terms).
%       Other options as taken by MAKEREGRESSOR including the important options
%        'noptions' and 'codeincr'.
% 
% See also MAKEREGRESSOR and FACT2REG

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
derivfun = [];
cliplimits = [-Inf Inf];
premult = 1;

while i <= length(varargin)
    
    switch lower(varargin{i})
        
%         case 'normconst' %multiply by this value before building polynomials.
%             normconst = varargin{i+1};
%             i = i+1;
        case 'label' %prepend this string to the regressor label
            label  = varargin{i+1};
            i = i+1; 
        case 'noptions' %number of options vector
            noptions  = varargin{i+1};
            i = i+1; 
        case 'premultiply' %multiply by this value before constructin the polynomial
            premult  = varargin{i+1};
            i = i+1; 
        case 'postmultiply' %after the regressor matrix is computed multiply each column by this vector
            postmult  = varargin{i+1};
            i = i+1; 
        case 'subtractdc' %subtracts a constant term equivalent to decorrelating with dc over the unit invterval
            subtractdc = varargin{i+1};
            i = i+1; 
        case 'codeincr' 
           codeincr = varargin{i+1};
            i = i+1; 
        case 'trig' %Polynomials are trigonometric, that is, sinusoids.
            trig = true;
        case 'dc' %include the "DC" zeroth order term (default is to start at 1st order)
            minord  = 0;
            subtractdc = false;
        case 'center' %Subtract mean value for each polynomial term except DC (if it's included)
            center  = true;
        case 'clip' %for values of the input outside this range, the output is  zero
            cliplimits  = varargin{i+1};
            i = i+1; 
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

dim = size(x,2);

intxnindx = zeros(1,dim*(1+dim)./2);


if ~trig    
    if dim == 1
        polyfunsdc = @(x,subdc) (repmat(x,1,length(poly)).^repmat(poly,size(x,1),1) - repmat(subdc,length(x),1))...
                        .*repmat(x>cliplimits(1),1,length(poly)).*repmat(x<=cliplimits(2),1,length(poly)); 
        polyfun = @(x) polyfunsdc(x*premult,subdc); 
        ncol = length(poly);
        
         Dmat = diag(1:polyord,1); %Derivative Matrix 
         if minord == 0
            derivfun = @(X,n) polyfunsdc(X*premult,subdc*(n==0))*Dmat^n;  %Function which returns the nth order polynomial derivative
        elseif minord == 1
            derivfun = @(X,n) ([(X*premult>cliplimits(1)).*(X*premult<=cliplimits(2)),polyfunsdc(X*premult,subdc*(n==0))]*Dmat^n)*diag(ones(1,polyord),-1)*eye(polyord+1,polyord); 
        end

    else
        ind = 1;
        if isscalar(polyord)
            polyord = repmat(polyord,1,dim);
        end
         
        ind = 1;
        ind2 = 1;
        startind = 0;
        endind = 0;
        Rintxn = [];
        argind = [];
        
        for i = 1:dim
            for j = 1:polyord(i)
                R(ind) =  buildpolyreg(x(:,i)*premult,j,'polyvec','subtractdc',false,varargin{:});
                R(ind).levmat(:) = j;
                R(ind).info.intxnord = j;
                R(ind).label = sprintf('x%i^%i',i,j);
                R(ind).code = ind;
                R(ind).codevec(1:R(ind).Npar) = ind;
                R(ind).factmat(:) = i;
                ind = ind+1;
            end
            if i > 1
                r2 = [];
                if ~isempty(Rintxn)
                    r2 = interaction(Rintxn,R(startind:end),'intxnord',polyord(i),'maxord',polyord(i));
                end
                Rintxn = [Rintxn, interaction(R(1:startind-1),R(startind:end),'intxnord',polyord(i),'maxord',polyord(i)),r2];
                
            end
                    
            startind = ind;

        end
        
        if minord == 0
            R0 = makeregressor(ones(size(x,1),1),'noptions',noptions);
            R0.factmat = 1; R0.levmat = 1;
            R0.function = @(X) ones(size(X));
        else
            R0 = [];
        end
        R = [R0,R,Rintxn];
        funs = {R.function};
        npars = [R.Npar];
        
        polygen = '@(x) [ ';
        for i = 1:length(R)
            indx = R(i).factmat(:,1);
            indx = indx(indx~=0);
            polygen = [polygen,sprintf('funs{%i}(',i),sprintf('x(:,%i),',indx),sprintf('),')];
        end

        polygen = regexprep(polygen,',)',')');
        polygen = [polygen,']'];
        polyfun = str2func(polygen);
        
        R = pool(R,'codeincr',codeincr,'label',label);
        
        %%% Creating matrix to compute derivatives...note this assumes that the
        %%% complete series of polynomials (1...p) are included...
        fmat = R.factmat';
        lmat = R.levmat';
        for i = 1:size(fmat,1)
            ordmat(i,fmat(i,fmat(i,:)~=0)) = lmat(i,fmat(i,:)~=0);
        end
        ordmat(end+1,:) = 0;
        
        for i = 1:dim
          dmat = zeros(R.Npar+1);
            dord = ordmat;
            dord(:,i) = dord(:,i)-1;
%             dord(dord<0) = 0;
            [ism,ismindx] = ismember(dord,ordmat,'rows');
            dmat(ismindx(ism),ordmat(:,i)~=0) = diag(ordmat(ordmat(:,i)~=0,i));
            Dmats{i} = dmat; %Dmats{k}^n*R.value computes the nth derivative wrt the kth variable 
        end
        
        
        
        if  subtractdc
            
             subdc = prod(1./(1+R.levmat),1); 
            R.function = @(x)polyfun(x*premult) - repmat(subdc,size(polyfun(x*premult),1),1);
            R.deriv = @(x,n) getderivx(polyfun(x*premult),n,Dmats) - (sum(n)==0)*repmat(subdc,size(polyfun(x*premult),1),1);
            
        else       
            R.function = polyfun;
            R.deriv = @(x,n) getderivx(polyfun(x*premult),n,Dmats);
            
            if minord == 0
                R.levmat(:,1) = 0;
            end
        end
%         R.factmat = R.code*ones(1,size(R.factmat,2));
%         R.levmat = 1:size(R.factmat,2);
        
        return
    end
        

else
     mkfarg = {};
    if ~subtractdc 
        mkfarg = {'keepdc'};
    end
    [makefourier_output,fr,deriv] = makefourier(polyord*ones(1, dim),mkfarg{:});     %trig polynomials may be in an arbitrary number of dimension based on length of polyord
    polyfun = @(x)makefourier_output(x*premult,1)';    
    derivfun = @(x,ns)deriv(x*premult,ns,1)';    
%     ncol = 2*length(poly)-sum(poly==0);
    ncol = length(polyfun(x(1,:)));
    
end    

if center
    
    regfun = @(x) polyfun(x) - ones(size(x,1),1)*mean(polyfun(x),1)*diag((minord:polyord)>0);
%     regfun = @(x) polyfun(x) - ones(size(x,1),1)*mean(polyfun(x),1)*diag((minord:polyord)>0);
else
    regfun = polyfun;
end

R = makeregressor(regfun(x).*repmat(postmult,1,ncol),'label',label,'normconst',normconst,'noptions',noptions,'codeincr',codeincr);

R.function = @(X) makeme(X,regfun,polyord);
if ~isempty(derivfun)
    R.deriv = @(X,n) makeme(X,derivfun,polyord,1,n);
end
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
function R = makeme(X,regfun,polyord,defaultval,varargin)

if size(X,2) < length(polyord)   
    error('Input matrix must have %i columns',length(polyord))
end

if nargin >3
    if isempty(varargin) || isempty(varargin{1})
        varargin{1} = defaultval;
    end
end
R = regfun(X,varargin{:});


%-----------------------
function DMX = getderivx(X,n,Dmats)

DM = Dmats{1}^n(1);

for k = 2:length(Dmats)
    DM = DM*Dmats{k}^n(k);
end
DM(:,end) =[];
X(:,end+1) = 1;

DMX = X*DM;
