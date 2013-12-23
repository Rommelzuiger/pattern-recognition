%SIMDISTM  Distance matrix between two data sets found from similarities
%
%   D = SIMDISTM (A,B)
%
% INPUT
%   A   NxK matrix or dataset
%   B   MxK matrix or dataset (optional; default: B=A)
%
% OUTPUT
%   D   NxM dissimilarity matrix or dataset; D in [0,1]
%
% DESCRIPTION
% Computes the distance matrix D between two sets of vectors: A and B.
% Distances between vectors X and Y are computed based on the similarity
% formula:
%     SIM(X,Y) = 2 * (X'Y) / (||X||^2 + ||Y||^2)
%     D(X,Y)   = (1 - SIM(X,Y))
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% DEFAULT
%   B = A
%
% SEE ALSO
% JACSIMDISTM, CORRDISTM, COSDISTM, LPDISTM, EUDISTM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function D = simdistm(A,B)

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


aa = sum(a.*a,2);
bb = sum(b.*b,2)';
D  = (2 * (a * b')) ./ (aa(:,ones(rb,1)) + bb(ones(ra,1),:));
D  = 1 - D;

% Check numerical inaccuracy
D (find (D < eps)) = 0;   % Make sure that distances are nonnegative
if bisa,
  D = 0.5*(D+D');         % Make sure that distances are symmetric for D(A,A)
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
