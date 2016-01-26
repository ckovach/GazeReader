function R = cov2corr(S)

%One-liner to convert variance-covariance matrix to correlation

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

R = diag(diag(S).^-.5)*S*diag(diag(S).^-.5);

