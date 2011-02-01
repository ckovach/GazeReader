function R = cov2corr(S)

%One-liner to convert variance-covariance matrix to correlation

R = diag(diag(S).^-.5)*S*diag(diag(S).^-.5);

