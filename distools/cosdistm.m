%COSDISTM  Distance Matrix based on Inner Products
%
%   D = COSDISTM(A,B)
%       OR
%   D = COSDISTM(A)
%
% INPUT
%   A 	NxK Matrix or dataset
%   B   MxK Matrix or dataset (optional; default: B=A)
%
% OUTPUT
%   D   NxM Dissimilarity matrix or dataset; D in [0,1]
%
% DESCRIPTION
% Computes a distance matrix D between two sets of vectors, A and B.
% Distances between vectors X and Y are derived based on their inner products
% (and their relations to the cosinus of the angle between them) as:
%     D(X,Y) = (1 - X'*Y/(||X||*||Y||))
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% DEFAULT
%   B = A
%
% SEE ALSO
% SIMDISTM, JACSIMDISTM, CORRDISTM, LPDISTM, EUDISTM
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function D = cosdistm (A,B)
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

if ~bisa,
  ab = a * b';
  a2 = sum(a.*a,2);
  b2 = sum(b.*b,2)';
  D  = ab ./ sqrt(a2(:,ones(rb,1)) .* b2(ones(ra,1),:));
else
  aa = a * a';
  a2 = diag(aa);
  a3 = diag(aa)';
  D  = aa ./ sqrt(a2(:,ones(ra,1)) .* a3(ones(ra,1),:));
end
D = (1 - D);
D(find (D < eps)) = 0;


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
