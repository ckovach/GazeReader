
function [lodds, loddserr,loddserrmat] = logoddserror(phis1,theta,S,phis2)

%  [lodds, loddserr, loddserrmat] = logoddserror(x,theta,S)
% Given the input regressor and asymptotic error covariance matrix, S, for a collection of bins,
% this function returns the mle estimate and estimation error for the log odds ratio over
% the bins, relative to uniform distribution, given input x. x is an MxK matrix 
% with regressor values for each of M bins. theta is a Kx1 vector of parameter values. 
% 
% [lpr, lprerr,lprerrmat] = logoddserror(x1,theta,S,x2)
% 
% Provides maximum likelihood and error estimate for the log odds ratio of
% outcome probabilities for ph1 over phi2. Note that x1 and x2 are both
% MxK matrices, where K is the value of each input regressor for each of M
% multinomial outcomes.
%
% See also ModelFit, mnlfit
%

% C. Kovach 2011



nbin = size(phis1,1);
npar = size(phis1,2);

rhos1 = phis1*theta;
pr1 = exp(rhos1)./nansum(exp(rhos1(:)));

phis1(pr1==0 | isnan(pr1),:) = 0; %This avoids possible nan if phis1 == +/-Inf
rhos1(pr1==0 | isnan(pr1),:) = -1./eps; 
pr1(isnan(pr1)) = 0; %treat nans as excluded bins
% exph1 = pr1'*phis1; %expected phi

if nargin > 3
    rhos2 = phis2*theta;
    pr2 = exp(rhos2)./nansum(exp(rhos2(:)));
    pr2(isnan(pr2)) = 0; %treat nans as excluded bins
    phis2(pr2==0 ,:) = 0; %This avoids possible nan if phis1 == +/-Inf
    rhos2(pr2==0 ,:) = -1./eps; 
    exph2 = pr2'*phis2; %expected phi
else
%     exph2 = 0;
    rhos2 = zeros(size(rhos1));
    phis2 = zeros(size(phis1));
end

prnot1 = zeros(nbin,nbin);
prnot2 = zeros(nbin,nbin);
for i = 1:length(rhos1)
    npr = zeros(1,nbin);
    npr([1:i-1,i+1:end]) = pr1([1:i-1,i+1:end]);
    prnot1(i,:) = npr./sum(npr);
    if nargin > 3
        npr2 = pr2([1:i-1,i+1:end]);
        prnot2(i,:) = npr2./sum(npr2);
    end
end


exph1not = prnot1*phis1; %Expectation of regressor over bins excluding the diag.
exph2not = prnot2*phis2;


C = phis1 - phis2 - (exph1not - exph2not);

lodds = rhos1 - rhos2 - (log(sum(exp(rhos1))-exp(rhos1))-log(sum(exp(rhos2))-exp(rhos2)));


Q = S*C';
loddserr = sqrt(sum(C'.*Q));
% lperr = sqrt(diag(lperrmat));

if nargout > 2
    loddserrmat = C*S*C';
end
