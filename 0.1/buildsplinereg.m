
function [R,Cmat] = buildsplinereg(X,knots,varargin)

% [R,Cmat] = buildsplinereg(X,knots,'splineord',splineord,'homogeneity',homogeneity)
% This function builds 1d polynomial regressors (R) and constraint matrix (C) for
% polynomial spline fitting using MODELFIT and MNLFIT. Knots are at the
% specified locations and the order of each splines is given by the
% 'splineord' argument (default is 3 for all). The parameter 'homogeneity' specifies
% the degree of continuity: degree n implies that the function is n-1 times
% differentiable at the knots. Default homogeneity is 2.
%
%SEE BUILDPOLYREG, MODELFIT, MNLFIT

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C Kovach 2010


homogeneity = 2;
nsplines = length(knots)-1;
return_pooled = true;
splineord = 3*ones(1,nsplines);    
cdi = 0;

i = 1;
bpargs = {};
while i <= length(varargin)
    
    switch lower(varargin{i})
        case 'homogeneity'                    % The parameter homogeneity specifies the degree of continuity: 
            homogeneity = varargin{i+1};     % degree n implies that the function is n times differentiable at the
            i = i+1;                            % knots.
        case 'splineord'                    
            splineord = varargin{i+1};    
            i = i+1;                            
        case 'splineord'                    
            splineord = varargin{i+1};    
            i = i+1;                            
        case 'dc'                    
            bpargs = [bpargs,{'dc'}] ;    
        case 'postmultiply'                    
            bpargs = [bpargs,{'postmultiply',varargin{i+1}}] ;    
            i = i+1;                            
        case 'codeincr'                    
            cdi = varargin{i+1};    
            i = i+1;                            
        otherwise
            bpargs = [bpargs,varargin{i},varargin{i+1}] ;    
            i = i+1;                            
    end
    i = i+1;
    
end



% sprange = [-Inf, knots,Inf];
clips = knots(cat(1,1:nsplines,2:nsplines+1))';

% subtr = clips(1,1);
% norm = 1;
    
%make polynomial regressors
for i = 1:nsplines   
    
    subtr(i) = clips(i,1);
    norm(i) = diff(clips(i,:));
    
    R(i) = buildpolyreg((X-subtr(i))./norm(i),splineord(i),'clip',[0 1],'codeincr',cdi,bpargs{:});
    
%         subtr(i+1) = clips(i+1,1);
%     if i < nsplines
%         norm(i+1) = diff(clips(i+1,1:2));
%     else
%         norm(i+1) = 1;
%     end
    
    cdi = cdi+1;
end

%make constraints
constraints = {};
totpars = sum([R.Npar]);
dt = 1e-9;
for i = 1:length(R)-1
    constraints{i} = zeros(homogeneity,totpars);
    D1 = R(i).function(1);
    D2 = R(i+1).function(0+eps);
%     D1 = R(i).function((knots(i)-dt);
%     D1 = R(i+1).function(knots(i));
    for h = 1:homogeneity-1
        D1(h+1,:) = R(i).deriv(1,h);
        D2(h+1,:) = R(i+1).deriv(0+eps,h);
%         D1(h+1,:) = R(i).deriv(knots(i)-dt,h);
%         D2(h+1,:) = R(i+1).deriv(knots(i),h);
    end
    
    constraints{i}(:,sum([R(1:i-1).Npar])+(1:R(i).Npar)) = D1;
    constraints{i}(:,sum([R(1:i).Npar])+(1:R(i+1).Npar)) = -D2;
end

Cmat = cat(1,constraints{:});

if return_pooled
    RP =  pool(R);

    RP.function = @(X) makefcn(X,RP.function,subtr,norm);
    RP.deriv = @(X,n) makefcn(X,RP.deriv,subtr,norm,n);
    RP.code = R(1).code;
    R =RP;
    
end
%%%

function RX = makefcn(X,fun,subtr,norm,varargin)

CX = {};

for i=  1:nargin(fun)./(length(varargin)+1)
    CX = cat(2,CX,{(X-subtr(i))./norm(i)},varargin);
end
    
%     CX(1:1+length(varargin):end) = repmat({(X-subtr)./norm},1,nargin(fun)./(length(varargin)+1));
%     for i = 1:length(varargin);
%         CX(i+1:1+length(varargin):end) = repmat(varargin(i),1,nargin(fun)./(length(varargin)+1));
%     end
    RX = fun(CX{:});
