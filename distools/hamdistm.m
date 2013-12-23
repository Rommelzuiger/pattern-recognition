%HAMDISTM Hamming Distance Matrix between Binary Vectors
%
%   D = HAMDISTM(A,B)
%       OR
%   D = HAMDISTM(A)
%
% INPUT
%   A   NxK Binary matrix or dataset
%   B   MxK Binary matrix or dataset
%
% OUTPUT
%   D   NxM Dissimilarity matrix or dataset
%
% DESCRIPTION
% Hamming distance between sets of binary vectors.
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function d = hamdistm(A,B)
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

if any(a~=0) | any(a~=1) | any(b~=0) | any(b~=1),
  error('Data should be binary.');
end

D = zeros(ra,rb);
for i=1:rb
  D(:,i) = sum((repmat(b(i,:),ra,1) ~= a),2);
end

% Set object labels and feature labels
if xor(isda, isdb),
  prwarning(1,'One matrix is a dataset and the other not. ')
end
if isda
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
