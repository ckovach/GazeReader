
function [rfunc,theta] = peak2gmm(peakdata, ij)

% [R,th0] = peak2gmm(peakdata,xx)
%
% This function takes a peakdata structure and returns 
% a function handle that generates regressors for a Gaussian mixture model, and
% returns the parameter to give corresponding peak locations and variances. 
%
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2010
%


if nargin < 2
    ij = false;
end

if ij == true
    posindx = [2 1];
else
    posindx = [1 2];
end


xx = [0 0];

npk = length(peakdata);

[pksrt,srti] = sort([peakdata.value],'descend');
peakdata = peakdata(srti);


rpoly = buildpolyreg(xx,2,'subtractdc',false);
rpolydc = buildpolyreg(xx,2,'dc','subtractdc',false);



funcs = {rpoly.function,rpolydc(ones(1,npk-1)).function}; 
    
rfunc = @(x,varargin) makereg(funcs,x,varargin{:});

eigv = peakdata(1).eigv(posindx,:);
eigs = peakdata(1).eigs(posindx);

H =eigv'*diag(eigs)*eigv; 

m = peakdata(1).pos(posindx)';

subth = cat(1,2*H*m,diag(-H),-2*H(1,2));
% subth = cat(1,H*m,diag(-H)/2,-H(1,2));

glwgt1 = peakdata(1).value - m'*H*m; 
theta{1} = subth([1 3 2 4 5]);



for i = 2:npk
    
    m = peakdata(i).pos(posindx)';

    eigv = peakdata(i).eigv(posindx,:);
    eigs = peakdata(i).eigs(posindx);
    H = eigv'*diag(eigs)*eigv; 
        
    glwgt = peakdata(i).value - m'*H*m; 
    
    
    subth = cat(1,2*H*m,diag(-H),-2*H(1,2));
%     subth = cat(1,H*m,diag(-H)/2,-H(1,2));
    
    theta{i} = cat(1,glwgt-glwgt1, subth([1 3 2 4 5]));
    
end

    

%%%%%%%%%%%%%%%5

function RX = makereg(funcs,xx,cdi,varargin)

%   RXs = cellfun(@(F) F(xx),funcs,'uniformoutput',false);
    
%   npars = cellfun(@(A) size(A,2),RXs);
  npk = length(funcs);
  
  for i = 1:npk
      v = zeros(npk,1);
      v(i) = 1;
      RX(i) = makeregressor(kron(funcs{i}(xx),v),'codeincr',cdi+i-1,'label',sprintf('gmm %i',i),varargin{:});
      RX(i).function = funcs{i};
  end
    