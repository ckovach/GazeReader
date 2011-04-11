
function [X,disc] = multitrim(X,trim,exvar);

% [X,Disc] = multitrim(X);
%Removes outliers from X using Multivariate iterative trimming of extreme data points. See Rencher p. 29.

if nargin < 2
    trim = .2;
end

if nargin < 3
  exvar = .95;
end

X = squeeze(X);



Xc = X;


mx = mean(Xc,2);
n = size(Xc,2);
%MX = mx(:,ones(size(Y,2),1));

Xc = Xc-mx*ones(1,size(Xc,2));

if size(X,1) > size(X,2) 
    S = (Xc'*Xc)/(n-1);   %Because S is singular, reduce dimensionality with svd.

    [U,D] = svd(S);

    d = diag(D);

    U = Xc*U*diag((d).^-.5);        %Converts eigenvectors of X'X to nonzero eigenvectors of XX'

    pckeep = min(find(cumsum(d/sum(d))>=exvar));

    r = rank(S);
    if pckeep > r | isempty(pckeep);
        pckeep = r;
    end

    d = d(1:pckeep);
    U = U(:,1:pckeep);

    Xcu = U'*Xc;
else
    Xcu = Xc;
end

zold = inf;
Disc = [];

Xu = Xcu;  
while 1

 Su = Xcu*Xcu';
 R = diag(diag(Su).^-.5)*Su*diag(diag(Su).^-.5);
 rs = R(find(triu(R+1e-15,1)));
 
 znew = .5*log((1+rs)./(1-rs));  % Fisher's z transform
  
 if mean(abs(znew-zold)) < .001
     break
 end
 
 zold = znew;
% T = diag(Xcu'*Su^-1*Xcu); 
  T = diag(Xu'*Su^-1*Xu); 

 [tsort,tind] = sort(T);
 disc = tind(end - ceil(trim*length(tind)) + 1:end);
 %Disc = cat(1,Disc,disc);
 %Xcu(:,disc) = [];
 Xcu = Xu; 
 Xcu(:,disc) = [];
 Xcu = Xcu - mean(Xcu,2)*ones(1,size(Xcu,2));
 if rank(Su) < max(size(Su));
     warning('Matrix is singular.')
     break
 end
end

X(:,disc) = [];

