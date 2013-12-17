function B = trainRidge(xtr,ytr,L)

%--------------------------------------------------------
% B = trainRidge(xtr,ytr,L)
%
% Trains a ridge regression model on the data in 
% (Xtrain,ytrain) with the penalty parameter L
%
% INPUT:
%  xtr: normalized data matrix
%  ytr: response variable
%    L: the penalty parameter, lambda
%
% OUTPUT:
%    B: the coefficents obtained with ridge regression
%--------------------------------------------------------

beta      = inv(xtr'*xtr + diag(ones(size(xtr,2),1),0).*L)*xtr'*(ytr-mean(ytr));
beta0     = mean(ytr);
B         = [beta0;beta];
return
