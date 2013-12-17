function B = trainLASSO(xtr,ytr,L)

%--------------------------------------------------------
% B = trainLASSO(xtr,ytr,L)
%
% Trains a LASSO regression model on the data in 
% (Xtrain,ytrain) with the penalty parameter L
%
% INPUT:
%  xtr: normalized data matrix
%  ytr: response variable
%    L: the penalty parameter, lambda
%
% OUTPUT:
%    B: the coefficents obtained with LASSO regression
%--------------------------------------------------------

[beta,msr] = arrfit(xtr,(ytr-mean(ytr)),L);
beta0      = mean(ytr);
B          = [beta0;beta];
return
