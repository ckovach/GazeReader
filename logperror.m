
function [lp, lperr,lperrmat,C] = logperror(phis1,theta,S,phis2)

%  [lp, lperr, lperrmat] = logperror(x,theta,S)
%
% Given the input regressor and asymptotic error covariance matrix, S, for a collection of bins,
% this function returns the mle estimate and estimation error for the log probability over
% the bins, given input x. x is an MxK matrix with regressor values for each
% of M bins. theta is a Kx1 vector of parameter values. 
% 
% [lpr, lprerr,lprerrmat] = logperror(x1,theta,S,x2)
% 
% Provides maximum likelihood and error estimate for the log ratio of
% outcome probabilities for ph1 over phi2. Note that x1 and x2 are both
% MxK matrices, where K is the value of each input regressor for each of M
% multinomial outcomes.
%
% [lpr, lprerr,lprerrmat,C] = logperror(x1,theta,S,x2)
% 
%  Also returns the contrast matrix which gives the first derivative of
%  lpr with respect to theta. Estimation error can be approximated as C*S*C', where C
%  is the covariance matrix for parameter estimates.
%
% See also ModelFit, mnlfit, logoddserror
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2010



nbin = size(phis1,1);
npar = size(phis1,2);

rhos1 = phis1*theta;
pr1 = exp(rhos1)./nansum(exp(rhos1(:)));
phis1(pr1==0 | isnan(pr1),:) = 0; %This avoids possible nan if phis1 == +/-Inf
rhos1(pr1==0 | isnan(pr1),:) = -1./eps; 
pr1(isnan(pr1)) = 0; %treat nans as excluded bins
exph1 = pr1'*phis1; %expected phi

if nargin > 3
    rhos2 = phis2*theta;
    pr2 = exp(rhos2)./nansum(exp(rhos2(:)));
    pr2(isnan(pr2)) = 0; %treat nans as excluded bins
    phis2(pr2==0 ,:) = 0; %This avoids possible nan if phis1 == +/-Inf
    rhos2(pr2==0 ,:) = -1./eps; 
    exph2 = pr2'*phis2; %expected phi
else
    exph2 = 0;
    rhos2 = 0;
    phis2 = 0;
end

C = phis1 - phis2 -repmat((exph1 - exph2),nbin,1);

lp = rhos1 - rhos2 - (log(sum(exp(rhos1)))-log(sum(exp(rhos2))));

Q = S*C';
lperr = sqrt(sum(C'.*Q));
% lperr = sqrt(diag(lperrmat));

if nargout > 2
    lperrmat = C*S*C';
end
