function Y = my_pdist(X,Y,s,t)
%PDIST Pairwise distance between observations.
%   Y = PDIST(X,METRIC) returns a vector which contains all the
%   distances between each pair of observations in X computed using
%   the given METRIC.  X is a M by N matrix, treated as M observations
%   of N variables. Since there are M*(M-1)/2 pairs of observations in
%   X, the size of Y is M*(M-1)/2 by 1.  The default metric is
%   'EUCLID'.  The available metrics are:
%
%      'euclid'    --- Euclidean metric
%      'seuclid'   --- Standardized Euclid metric
%      'cityblock' --- City Block metric
%      'mahal'     --- Mahalanobis metric
%      'minkowski' --- Minkowski metric
%
%   Y = PDIST(X, 'minkowski', p) specifies the exponents in the
%   'Minkowski' computation. When p is not given, p is set to 2.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,M),
%   (2,3),...(2,M),.....(M-1,M)).  i.e. the upper right triangle of
%   the M by M square matrix. To get the distance between observation
%   i and observation j, either use the formula Y((i-1)*(M-i/2)+j-i)
%   or use the helper function Z = SQUAREFORM(Y), which will return a
%   M by M symmetric square matrix, with Z(i,j) equaling the distance
%   between observation i and observation j.
%
%   See also SQUAREFORM, LINKAGE

%   ZP You, 3-10-98
%   Copyright 1993-2000 The MathWorks, Inc.
%   $Revision: 1.1 $

if nargin >= 3
   if length(s) < 2
      error('Unrecognized metric');
   else 
      s = lower(s(1:2));
   end
else
   s = 'eu';
end

if s == 'mi' % Minkowski distance need a third argument
   if nargin < 4
      t = 2; 
   elseif t <= 0
      error('The fourth argument has to be positive.');
   end
end

[mx, nx] = size(X);
[my, ny] = size(Y);

%if (mx < 2)
%   error('The first argument has to be a numerical matrix with at least two rows');
%end

I = (1:mx)'*ones(1,my); I = reshape(I,mx*my,1);
J = (1:my)'*ones(1,mx); J = reshape(J',mx*my,1);

Y = (X(I,:)-Y(J,:))';
I = []; J = []; p = [];  % no need for I J p any more.

switch s
case 'eu' % Euclidean
   Y = sum(Y.^2,1);
   Y = sqrt(Y);
case 'se' % Standadized Euclidean
   D = diag(1./var(X));
   Y = sum(D*(Y.^2),1);
   Y = sqrt(Y);
case 'ci' % City Block
   Y = sum(abs(Y),1);
case 'ma' % Mahalanobis
   v = inv(cov(X));
   Y = sqrt(sum((v*Y).*Y,1));
case 'mi' % Minkowski
   Y = sum(abs(Y).^t,1).^(1/t);
otherwise
   error('no such method.');
end
