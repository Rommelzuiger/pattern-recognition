%CHECKTR Check whether Square Distance Matrix Obeys the Triangle Inequality
%
%   [M,C] = CHECKTR(D)
%
% INPUT
%   D   NxN symmetric distance matrix or dataset
%
% OUTPUT
%   M   Total number of inequalities disobeying the triangle inequality
%   C   Constant determined such that when added to all off-diagonal
%       values of D, makes D metric
%
% DESCRIPTION
% Checks whether a symmetric distance matrix D obeys the triangle inequality.
% If D is asymmetric, it is made symmetric by averaging. M is the total number
% of disobeyed inequalities. Hence, for M = 0, the triangle inequality holds
% for the distance matrix D. Additionally, a constant C is sought, which when
% added to the off-diagonal elements of D will make D obey the triangle inequality:
%    C = max_{i,j,k} {abs(D_{ij} + D_{ik}  - D_{jk})}
% If C is zero, then D fulfills the condition.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [no,c] = checktr(d)
d = +d;
[n,m] = size(d);

prec = 1e-12;
if ~issym(d,prec)
  prwarning(1,'Distance matrix should be symmetric. It is made so by averaging.')
end
d = 0.5*(d+d');
d(1:n+1:end) = 0;

c = 0;
ismetric = 1;

no = 0;         % Number of disobeyed triangle inequalities
for j=1:n-1
  r1 = d(:,j);
  r2 = d(j,j+1:n);
  % M checks the triangle ineqaulities; all elements of M are
  % positive or zero if this holds.
  M  = repmat(r1,1,n-j) + d(:,j+1:n) - repmat(r2,n,1);

  % The elemnets of M should ideally be compared to 0. However,
  % due to numerical inaccuracies, a small negative number should
  % be used. Experiments indicate that -1e-13 works well.
  % Otherwise, one gets wrong results.
  tol = 1e-13;
  mm  = sum(M(:) < -tol);

  no = no + mm;
  ismetric = ismetric & (mm == 0);
  if ~ismetric,
    cc = max(max(M));
    if cc > c,
      c = cc;
    end
  end
end
return;
