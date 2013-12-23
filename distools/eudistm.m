%EUDISTM Euclidean Distance Matrix 
%
%   D = EUDISTM(A,B)
%       OR
%   D = EUDISTM(A)
%
% INPUT
%   A   NxK Matrix or dataset
%   B   MxK Matrix or dataset (optional; default: B = A)
%
% OUTPUT
%   D   NxM Euclidean distance dataset or matrix
%
% DESCRIPTION
% Computation of the Euclidean distance matrix D between two sets of vectors
% A and B. If A and B are datasets, then D is a dataset as well with the labels
% defined by the labels of A and the feature labels defined by the labels of B.
% If A is not a dataset, but a matrix of doubles then D is also a matrix of
% doubles.
%
% NOTE
% EUDISTM(A,B) is equivalent to SQRT(DISTM(A,B)).
%
% DEFAULT
%   B = A
%
% SEE ALSO
% DATASETS, DISTM, PROXM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function D = eudistm(A,B)
bisa = nargin < 2;
if bisa,
  B = A;
end

isda = isdataset(A);
isdb = isdataset(B);
a = +A;
b = +B;

[ra,ca] = size(a);
[rb,cb] = size(b);
if ca ~= cb,
  error ('Matrices should have equal numbers of columns');
end


% The order of operations below is good for the accuracy.
D = ones(ra,1)*sum(b'.*b',1);
D = D + sum(a'.*a',1)'*ones(1,rb);
D = D - 2 .*(+a)*(+b)';

% Check for a numerical inaccuracy
D(find(D<eps)) = 0;

D = sqrt(D);

% Take care of symmetric distance matrix
if bisa & ra == rb,
  D = 0.5*(D + D');
  D([1:ra+1:ra^2]) = 0;
end

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
