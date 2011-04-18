
function peakdata = findpeaks(dfunc,parest,I, range, density,type,fig)


%  peakdata = findpeaks(dfunc,parest,I, range, density,type)
%
% Find peaks, troughs and saddle points over a given observation space.
% The bounds of the space are given by range, where the range for the 
% kth dimension is range(k,1) - range(k,2). Initially the space is sampled
% over a regular grid with density k.
%
% Specifying type returns vanishing points of the given type: 1 = peaks, 0
% = saddle points, -1 = troughs. 
%
% peakdata is a structure with information on each of the returned points.
%
% peakdata.eigs contains the greatest and least eigenvaues of the 2nd
% derivative tensor.
%
% peakdata.eigv contains the corresponding eigenvectors.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


% C Kovach 2010
% 

ndim = 2;

if nargin < 7 || isempty(fig)
    makeplot = false; 
else
    makeplot = true;
end


plotcyc = false; 

if nargin < 4 || isempty(range)
    range = repmat([0 1],ndim,1);
end

if nargin < 5 || isempty(density)
    density = 50*ones(1,ndim);
end
if nargin < 6 
	type = 1;
end

for i = 1:ndim
    ldx{i} = linspace(range(i,1),range(i,2),density(i));
end

xs = cell(1,ndim);

if ndim == 2
    ord = [2 1];
else
    ord = 1:ndim;
end

[xs{ord}] = ndgrid(ldx{:});
xx = cellfun(@(a) a(:), xs,'uniformoutput',false);
xx = [xx{:}]; 
npt = size(xx,1);

xx0 = xx;

stthresh = 0;


tol = 1e-6;

dxx = Inf;

thresh = 1; %discard point if increment is larger than this value
xxold = Inf(npt,ndim);

if plotcyc
    fgc = figure;
    plc = plot(xx0(:,1),xx(:,2),'.');
    axis([0 1 0 1])
    
end

cyclehist = 10; 
maxiter = 500;
iter = 0;
disc = false(npt,1);
while max(sqrt(sum(dxx.^2,2))) > tol && iter < maxiter
    Ds = zeros([ndim npt]);
    D2s = zeros(ndim,ndim*npt);
    for i = 1:ndim
        %First deriv
        nd = zeros(1,ndim); 
        nd(i) = 1;
        ddx = dfunc(xx,nd);
        Ds(i,:) = ddx*parest;

        %Second derivative
        for j = 1:ndim
            nd2 = nd; 
            nd2(j) = nd2(j) + 1;        
            ddx2 =  dfunc(xx,nd2);
            D2s(i,j+(0:size(xx,1)-1)*ndim) = ddx2*parest;
        end

    end

    D2sp = sparseblock(D2s,ndim);
    
    dxx = reshape(-D2sp\Ds(:),ndim,npt)';
    
    
   
    
    if plotcyc
        set(plc,'xdata',xx(:,1),'ydata',xx(:,2));
        drawnow;
    end
    
    xx = xx+dxx;
    
    disc = all( sqrt(sum((xxold-xx).^2,2)) < 1e-10 ,2) & sqrt(sum(dxx.^2,2)) > tol;  %To avoid unending cycles, all points which repeat a previous (other than most recent) position are discarded  
    
    disc = disc  |  any(xx<repmat(range(:,1)',npt,1),2)  |  any(xx>repmat(range(:,2)',npt,1),2);
    
    xx = xx(~disc,:);

    if mod(iter,cyclehist) == 0 || any(disc)
        xxold = xx;
    end
    
    iter = iter +1;
    
%     xxs(iter,:) = xx;
    
    
    npt = size(xx,1);
%     max(abs(dxx(:)))
end

if iter >= maxiter, warning('Maximum iterations reached\npeak identification may be inaccurate'), end

 Ds = zeros([ndim npt]);
for i = 1:ndim
    %First deriv
    nd = zeros(1,ndim); 
    nd(i) = 1;
    ddx = dfunc(xx,nd);
    Ds(i,:) = ddx*parest;

end

smDs = sqrt(sum(Ds.^2));

xx(smDs>tol,:) = []; %exclude points which didn't converge to tolerance (these are probably cyclic).

utol = 1e-3;
% utol = .01;
[uqx,uqi] = unique(round(xx./utol)*utol,'rows'); %These are the unique points at which the derivative vanishes: peaks, troughs and saddles
                                                 %Next we obtain second
                                                 %order statistics to decide which

[eigs,eigv, wald] = deriv2(dfunc,parest,I,uqx);

npks = size(eigs,2);


for i = 1:npks
    
    peakdata(i).pos = uqx(i,:);
    peakdata(i).value = dfunc(uqx(i,:),[0 0])*parest;
    peakdata(i).eigs = eigs(:,i);
    peakdata(i).eigv = squeeze(eigv(:,i,:));    
    peakdata(i).type = (eigs(1,i) > 0) - (eigs(2,i) < 0);
    peakdata(i).walds = wald(:,i);
    
    if peakdata(i).type ~=0 
        peakdata(i).typestat = wald( 3/2 - 1/2*peakdata(i).type,i);
    else
        peakdata(i).typestat = nan;
    end        
    
    
end

if ~isempty(type)
    peakdata = peakdata([peakdata.type] == type);
end

if makeplot && ndim == 2
    
    if nargin < 7
        fig = figure; 
    end
    rho = reshape(dfunc(xx0,[0 0])*parest,density);
    imagesc(range(1,:),range(2,:),exp(rho))
    hold on

    
    pl = plotpeaks(peakdata(abs([peakdata.typestat])>stthresh & [peakdata.type] ~=0),fig,[],'k');
end






                                                 

