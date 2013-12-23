%CORRTR Correct square distance matrix for the violation of the triangle inequality
%
%   DC = CORRTR(D)
%
% INPUT
%   D   NxN symmetric dissimilarity matrix or dataset
%
% OUTPUT
%   DC  Dissimilarity matrix/dataset with all entries corrected where the triangle inequality is violated
%

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org, Ela Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function d = corrtr(x)

d = +x;
[n,m] = size(d);
if n ~= m
	error('Dissimilarity matrix should be square')
end

d = min(d,d');
d(1:n+1:end) = 0;

for j=1:n-1 % correct the upper triangle, top down
  r1 = d(:,j);
  r2 = d(j,j+1:n);
  % M checks the triangle ineqaulities; all elements of M are
  % positive or zero if this holds.
  % M  = repmat(r1,1,n-j) - d(:,j+1:n) + repmat(r2,n,1);
  d(:,j+1:n)  = min(repmat(r1,1,n-j) + repmat(r2,n,1), d(:,j+1:n));
end
d = min(d,d'); % correct the bottom triangle as well

for j=n:-1:2 % correct the bottom triangle, bottom up
  r1 = d(:,j);
  r2 = d(j,1:j-1);
  % M checks the triangle ineqaulities; all elements of M are
  % positive or zero if this holds.
  % M  = repmat(r1,1,n-j) - d(:,j+1:n) + repmat(r2,n,1);
  d(:,1:j-1)  = min(repmat(r1,1,j-1) + repmat(r2,n,1), d(:,1:j-1));
end
d = min(d,d'); % correct the upper triangle as well

d = setdata(x,d);

return
