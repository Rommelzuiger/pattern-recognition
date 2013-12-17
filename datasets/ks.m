%KS Kernel smoothing
%
%    YTEST_HAT = KS(XTRAIN,YTRAIN,XTEST,H)
%
% Fit a kernel smoother on training data XTRAIN, YTRAIN and width
% H, and evaluate it on the test data XTEST.
function ytest_hat = ks(xtrain,ytrain,xtest,h)

for i=1:length(xtest),
   K = exp(-(xtrain-xtest(i)).^2/h.^2);
   ytest_hat(i,1) = sum(K .* ytrain) / sum(K);
end
