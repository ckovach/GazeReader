
function [PC,lc,AIC,BIC] = mlpca(X,varargin)

%Maximum likelihood pca with dimensionality chosen based on AIC or BIC
%See Tipping 1999

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


i = 1;

varargin{end+1} = 'finis';

criterion = 'bic';
center = false;
correl = false;

maxeigs = Inf;
maxdim = 300; %Maximum dimensionality
trim = [];
while i <= length(varargin)
    switch lower(varargin{i})
       case 'bic'
            criterion = 'bic';
       case 'aic'
            criterion = 'aic';
       case 'trim'
            trim  = varargin{i+1};
            i = i+1;
       case 'center'
            center= true;
       case 'correl'
            correl= true;
       case 'maxdim'
            maxdim = varargin{i+1};
            i = i+1;
       case 'finis'
         
      otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end 


if ~isempty(trim)
  [X,dsc] = multitrim(X,trim);
  factor(dsc,:) = [];
else
    trim = 0;
end

flip = 0;

if size(X,1) >= size(X,2) 
  S = X'*X;
  flip = 1;
else
  S = X*X';
end


if correl
    DS = diag(S).^.5;
    S = diag(DS.^-1)*S*diag(DS.^-1); %correlation matrix
end

if center && flip
    mp = mean(X,1);
    S = S - mp'*mp;         %Covariance matrix
elseif center
    mp = mean(X,2);
    S = S - mp*mp';
end

opts.issym = true;
opts.disp = 0;


totvar = trace(S);

K = length(S);
if K < maxdim
    neigs = K ;
    maxdim = K;
% else
%     neigs = maxeigs;
    
end

% [U,D] = eigs(S,neigs,'LA',opts);
[U,D] = svd(S);
l = diag(D);

sigml = (totvar - cumsum(l))'./(K-(1:length(l))); %ml estimate remaining variance

N = size(X,2);

%find lowest AIC or BIC

for i = 1:maxdim
    i
%     sg2(i) = 1./(K-i)*sum(l(i+1:end));
    
    W = U(:,1:i)*(diag((l(1:i)-sigml(i)).^.5));
    C = W*W'+sigml(i)*eye(K);
    
%     lc = [l(1:i);sigml(i)*ones(K-i,1)];
    lc = [l(1:i);zeros(K-i,1)]+sigml(i);
    
%     V = sum(X.^2,2);
    
    Cchol = chol(C);
    
    %log likelihood
%     a(i) = trace(Cchol\(Cchol'\S)); %b(i) = trace(Cchol^-1*(S/Cchol'^-1));
    LL = -N/2*(K*log(2*pi) + trace(Cchol\(Cchol'\S)) + sum(log(lc)));
    ll(i) = LL;
    %AIC

    npar = i*length(S) + 1 + i*(i-1)./2;
%     npar = 1+i*(i+1)/2+i;
    
%     AIC(i) = -2*LL + 2*i  + 2*i*(i+1)./(N-i-1);
    AIC(i) = -2*LL + 2*npar;
    
    %     AICc(i) = -2*LL + 2*npar  + 2*npar*(npar+1)./(N-npar-1);
    
    %BIC
%     BIC(i) = -2*LL + log(N)*i;
    BIC(i) = -2*LL + log(N)*npar;
   
     ICPCA(i) = sum(log(lc)) + i./N*log(N);
       ICPCA2(i) = sum(log(lc)) + npar./N*log(N);
    
    np(i) = npar;
end
    
if strcmp(criterion,'bic')
    [mn,npc] = min(BIC);
elseif strcmp(criterion,'aic')
    [mn,npc] = min(AIC);
end   

  
lc = [l(1:npc);sigml(npc)*ones(K-npc,1)];

if flip
    U = X*U(1:npc)*diag(l(1:npc).^-.5);        %Converts eigenvectors of X'X to nonzero eigenvectors of XX'
end

PC = U(:,1:npc);

    
