%EXPDISTM  Exponential-type of Distance Matrix
%
%     D = EXPDISTM (A,B,R)
%       or
%     D = EXPDISTM (A,B)
%       or
%     D = EXPDISTM (A,R)
%       or
%     D = EXPDISTM (A)
%
% INPUT
%   A   NxK Matrix or dataset
%   B   MxK Matrix or dataset (optional; default: B=A)
%   R   Parameter to scale the Gaussian function
%       (optional; default: sqrt(N)*(max(max(A))- min(min(A))))
%
% OUTPUT
%   D   NxM Dissimilarity matrix or dataset
%
% DESCRIPTION
% Computes the distance matrix D between two sets of vectors, A and B.
% Given the vectors X and Y, distances are computed as:
%   D(X,Y) = 1 - exp (-(X-Y)'(X-Y)/R^2)
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% SEE ALSO
% SIMDISTM, JACSIMDISTM, COSDISTM, CORRDISTM, LPDISTM, EUDISTM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [D,r] = expdistm(A,B,r)
[ra,ca] = size(A);

bisa = 0;
if nargin < 2,
  bisa = 1;
  r = (max(max(+A)) - min(min(+A)))*sqrt(ra);
  B = A;
  [rb,cb] = size(B);
else
  [rb,cb] = size(B);

  if nargin < 3,
    if max (rb,cb) == 1,
      bisa = 1;
      r = B;
      B = A;
      [rb,cb] = size(B);
    else
      r = (max(max(+A)) - min(min(+A)))*sqrt(ra);
    end
  end
end


if ~bisa & ca ~= cb,
  error ('The matrices should have the same number of columns.');
end

isda = isa(A,'dataset');
isdb = isa(B,'dataset');
a = +A;
b = +B;


%d = zeros(ra,rb);
%d = sum ((abs (repmat (permute(a,[1 3 2]), [1 rb 1]) - ...
%               repmat (permute(b,[3 1 2]), [ra 1 1]))).^2,3);

D = 1 - exp (-distm(a,b)/r^2);
D(find(D < eps)) = 0;


% Set object labels and feature labels
if xor(isda, isdb),
  prwarning(1,'One matrix is a dataset and the other not. ')
end
if isda,
  if isdb,
    D = setdata(A,D,getlab(B));
  else
    D = setdata(A,D);
  end
  D.name = 'Distance matrix';
  if ~isempty(A.name)
    D.name = [D.name ' for ' A.name];
  end
end
return
