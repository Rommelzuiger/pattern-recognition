%TRAINLINREG Linear regression
%
%   [beta0_hat,beta1_hat] = trainlinreg(xtrain,ytrain)
%
% Compute the least squares solution for finding the linear regressor
% through the data (xtrain,ytrain), where xtrain and ytrain are Nx1
% vectgors. The slope beta1 and offset beta0 are returned.
function [beta0_hat,beta1_hat] = trainlinreg(xtrain,ytrain)

if size(xtrain,1)==1
	xtrain = xtrain';
end
if size(ytrain,1)==1
	ytrain = ytrain';
end

X = [ones(size(xtrain,1),1) xtrain];
B = inv(X'*X)*X'*ytrain;
beta0_hat = B(1);
beta1_hat = B(2:end);

return
