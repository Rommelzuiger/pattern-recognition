function [mse,yte_hat] = testlinreg(B,xte,yte)

%--------------------------------------------------------
% [mse,yte_hat] = testlinreg(B,xte,yte)
%
% Tests a linear regression model with coefficients
% in B on the data in (xte,yte) and computes the
% mean squared error of the estimated outputs
% with respect to the true outputs.
%
% INPUT:
%    B: ridge regression coefficients
%  xte: normalized data matrix
%  yte: response variable
%
% OUTPUT:
%   mse: mean squared error of the estimated outputs
%        with respect to the true outputs.
%  yte_hat: estimated outputs (responses)
%--------------------------------------------------------

beta0   = B(1);
beta    = B(2:end);
yte_hat = beta0 + xte*beta;
mse     = mean((yte_hat - yte).^2);
return
