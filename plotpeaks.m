
function pl = plotpeaks(peakdata,fig,xyscale,varargin)
% function pl = plotpeaks(peakdata,fig,ij,varargin)

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% plotpeaks(peakdata,fig,ij,varargin)
%
%Plots a standard error ellipse at peak locations.
%
% If ij = true (default), then the x- and y- axes are reversed in the plot.
%
% Subsequent arguments are passed to PLOT.
%
%

% C. Kovach 2010



if nargin < 2 || isempty(fig)
    fig = figure;
end

ij = false;

figure(fig)

hold on,

if nargin < 3 || isempty(xyscale)
    xyscale = 1;   
end

if ij
    axind = [2 1];
else
    axind = [1 2];
end

for i = 1:length(peakdata)
        if length(peakdata(i).pos) ~=2, warning('This function plots peaks only in 2 dimensions'); continue, end
        if peakdata(i).type ~=0
            if peakdata(i).type == 1 && nargin < 4
                varargin = {'m'}; elseif peakdata(i).type == -1 && nargin < 4 varargin = {'c'};
            end
            
            eigv = peakdata(i).eigv(axind,:);
            eigs = diag(xyscale.^-2)*peakdata(i).eigs;
            
            T = eigv*diag(( peakdata(i).type*eigs).^-.5)*eigv';

            pl(i,:) = ellipse(peakdata(i).pos(axind)*diag(xyscale),T ,fig,varargin{:});
        else
            warning('Saddle points are not plotted')
        end
        
end