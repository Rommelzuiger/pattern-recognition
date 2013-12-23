%ISSYM Check whether a Matrix is Symmetric
%
%   OK = ISSYM(A,DELTA)
%
% INPUT
%   A      Dataset
%   DELTA  Parameter for the precision check (optional; default: 1e-13)
%
% OUTPUT
%   OK     1 if the matrix A is symmetric and error otherwise.
%
% DESCRIPTION
% A is considered as a symmetric matrix, when it is square and
% max(max(A-A')) is smaller than DELTA.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [ok,nn] = issym(A,delta)

if nargin < 2,
  prwarning(6,'The precision is not provided, set up to 1e-13.');
  delta = 1e-13;
end

A = +A;
[m,k] = size(A);
if m ~= k,
  error ('Matrix should be square.')
end
nn = max(max((A-A')));
ok = (nn < delta);

if ~ok & nargout == 0
  error('Symmetric matrix expected.')
end
	
return;
