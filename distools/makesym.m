%MAKESYM Make a matrix symmetric
%
%   [B,C] = MAKESYM(A)
%
% INPUT
%   A   Dataset or matrix
%
% OUTPUT
%   B   Symmetric dataset or matrix computed as (A+A')/2
%   C   Asymmetric remaining part, dataset or matrix computed as (A-A')/2
%
% DESCRIPTION
% B is a symmetric matrix obtained by averiging values of A.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% EWI Faculty, Delft University of Technology and
% School of Computer Science, University of Manchester


function [d1,d2] = makesym(d)
[m,k] = size(d);
if m ~= k,
  error ('Matrix should be square.')
end
d1 = 0.5 * (d + d');
if nargout == 2,
  d2 = 0.5 * (d - d');
end
return
