%SIGMOID Element-wise Sigmoid Tranformation of a Matrix
%
%   F = SIGMOID(A,S)
%
% INPUT
%   A   NxM matrix (dataset)
%   S   Scale parameter (optional, default: 1)
%
% OUTPUT
%   F   NxM matrix (dataset)
%
% DESCRIPTION
% Applies sigmoid transformation to the elements of A.
% The values of F are in (0,1] for positive values of A and
% in [-1,0) for negative values.
%
% DEFAULT
%   S = 1
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function F = sigmoid(A,s)

if nargin < 2,
  s = 1;
end

F = 2./(1+exp(-A/s)) - 1;
return;
