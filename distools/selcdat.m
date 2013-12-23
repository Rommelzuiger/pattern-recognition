%SELCDAT Select Class Subset from a Square Dissimilarity Dataset
%
%   [DN,J] = SELCDAT(D,C)
%
% INPUT
%   A   NxN Dissimilarity Dataset
%   C   Indices of classes
%
% OUTPUT
%   DN  Subset of the dataset D
%   J   Indices of the selected objects
%
% DESCRIPTION
% The classes listed in C (numerically) are extracted for the square
% dissimilarity matrix D by both, their rows (objects) as well as their
% columns (features).

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org, and
% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [D,J] = selcdat(D,n)
issquare(D);
if nargin < 2, return, end

J = [];
c = getsize(D,3);
if (any(n > c))
  error('Not that many classes')
end

for j=1:length(n)
  J = [J; findnlab(D,n(j))];
end
D = remclass(D(J,J));
return;
