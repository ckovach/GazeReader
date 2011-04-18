
function [eigs,eigv, wald,D2s] =deriv2(dfunc,parest,Imat,range, density)

% [eigs,eigv, wald,S] = deriv2(dfunc,parest,I, xx)
%
% Returns second  derivative of the log density at each row in xx in the
% form of eigenvalues and eigenvectors of the 2nd derivative tensor. If the Fisher information
% matrix, I, is non-empty, then the wald statistic for minimum and maximum eigenvectors is returned
% in 'wald'.
%
% dfunc is a function that computes the analytic derivative at the specified points. For a
% regressor, R, with defined derivative function, dfunc = R.deriv.
%
% [eigs,eigv, wald] = deriv2(dfunc,parest,I, range,density)
% Returns values at a regularly sampled grid, where dimension k is sampled
% from range(k,1) to range(k,2) at density(k).
% 
% see also FINDPEAKS
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2010



ndim = 2;
xx = [];
if nargin < 4 %range of sampling
    range = repmat([0 1],ndim,1);
    density = 100*ones(1,ndim);
end

if nargin == 4 || isempty(density) %density of sampling
%     density = 100*ones(1,ndim);
     xx = range;
     density = [];
end

if nargin < 3
 I = [];
end

if ~isempty(density)
    for i = 1:ndim
        ldx{i} = linspace(range(i,1),range(i,2),density(i));
    end
end

if isempty(xx)
    xs = cell(1,ndim);
    
    if ndim == 2
        ord = [2 1];
    else
        ord = [1 2];
    end
    
    [xs{ord}] = ndgrid(ldx{:});
    xx = cellfun(@(a) a(:), xs,'uniformoutput',false);
    xx = [xx{:}]; 
end



D2s = [];
for i = 1:ndim
    %First deriv
    nd = zeros(1,ndim); 
    nd(i) = 1;
    ddx = dfunc(xx,nd);
    Ds{i} = ddx*parest;
%     E = ddx*S;
%     Ws{i} = Ds{i}./sqrt(sum(E.*ddx,2)); %Wald statistic       
    
    %Second derivative
    for j = 1:ndim
        nd2 = nd; 
        nd2(j) = nd2(j) + 1;        
        ddx2 =  dfunc(xx,nd2);
        ddxs{i,j} = ddx2;
        
        D2s(i,j+(0:size(xx,1)-1)*ndim) = ddx2*parest;
    end
       
end

D2sp = sparseblock(D2s,ndim);

%Get eigenvalues and vectors....


npt = size(xx,1);

th0 = ones(ndim+1,1);
th0(1:ndim) = th0(1:ndim)./sqrt(ndim);


th0 = repmat(th0,1,npt);

lm = false(size(th0));
lm(ndim+1,:) = true;

lmq = false(ndim+1, numel(th0));
lmq(:,lm) = true;
% lmq = sparseblock(lmq,ndim+1)>0;
dgs = speye(size(D2sp))==1;
lmexpand = (kron(1:npt,ones(1,ndim))*ndim-1)';
lambda = spalloc(size(th0,1),size(th0,2),ndim*npt);

tol = 1e-6;

th = zeros(size(th0));    
%       dir = (2*mod(i,2)-1); 
      dir = 1; 
    QM = zeros(ndim+1, numel(th0));
    QM(1:ndim,~lm) = dir*2*D2s;
for i = 1:ndim  % finding the eigenvectors..
    
    th1 = th;
    th = th0;
    th(1:ndim,:) = th0(1:ndim,:) - repmat(sum(th1(1:ndim,:).*th0(1:ndim,:)),ndim,1).*th1(1:ndim,:);
    
    dth = Inf;

    while max(abs(dth(:))) > tol
        
        
        QM(1:ndim,lm) = dir*2*th(1:ndim,:);        
        QM(ndim+1,~lm) = dir*2*th(~lm);                 
        QM(ndim+1,lm) = 0;                 
        
        QMsp = sparseblock(QM,ndim+1);
        
        lambda(1:ndim,:) = th(ndim+ones(ndim,1),:);
        lambda(3,:) = 0;
        
        Q = QMsp*th(:);
%         Q = (2*mod(i,2)-1)*QMsp*th(:);
        Q(ndim+1:ndim+1:end) = dir*(sum(th(1:ndim,:).^2)-1);
%         Q(1:ndim,:) = 2*Q(1:ndim,:) + 2*lambda(1:ndim,:).*th(1:ndim,:);
%         Q(lm(:)) = sum(QM(:,lm).^2)-1;
        
%         QM(ndim+1,~lm(:)) = 2*th(~lm(:));
        
        DQ = QMsp + dir*2*diag(sparse(lambda(:)));
        
                
        dth = reshape(-DQ\Q,ndim+1,npt);
        
        th = th+dth;
    end
    
    ths(:,:,i) = th;
               
end

%Sort eigenvectors by  eigenvalue sign and magnitude
ls = squeeze(ths(ndim+1,:,:));
if npt > 1
    ls = ls';
end
[mx,mxi] = sort(ls,1);

[J,I] = meshgrid(1:size(ths,2),1:size(ths,1));
lind = sub2ind([size(th,1),size(th,2)],I,J);
linds = repmat(lind,[1 1 ndim]) + repmat(permute(mxi-1,[3 2 1])*numel(lind),[ndim+1 1 1]);
ths = ths(linds); %sorted

eigs = squeeze(ths(ndim+1,:,:));
if npt > 1
    eigs = eigs';
end
eigv = ths(1:ndim,:,:);



if ~isempty(Imat)
    S = Imat^-1;

    npar = length(Imat);
    dm = [1 ndim];
    for k = 1:2
        es = eigs(dm(k),:);
        ev = eigv(:,:,dm(k));

        econtrast = 0;  
        for i = 1:ndim
            for j = 1:ndim

                econtrast =  econtrast + ddxs{i,j}.*repmat((ev(i,:).*ev(j,:))',1,npar);
            end
        end
        
        err = sqrt(sum(econtrast.*(econtrast*S),2));
        wald(k,:) = (econtrast*parest)./err;
    end
end

% 
%     
