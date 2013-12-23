%CORRDISTM  Distance Matrix based on Correlations
%
%   D = CORRDISTM(A,B)
%       OR
%   D = CORRDISTM(A)
%
% INPUT
%   A   NxK Matrix or dataset
%   B   MxK Matrix or dataset (optional; default: B=A)
%
% OUTPUT
%   D   NxM dissimilarity matrix or dataset; D in [0,1]
%
% DESCRIPTION
% Computes the distance matrix D between two sets of vectors, A and B.
% Distances between vectors X and Y are computed based on their correlation as:
%     CORR(X,Y) = COV(X,Y) / sqrt(Var(X) * VAR(Y))
%     D(X,Y)    = (1 - CORR)/2
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% DEFAULT
%   B = A
%
% SEE ALSO
% SIMDISTM, JACSIMDISTM, COSDISTM, LPDISTM, EUDISTM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function D = corrdistm (A,B)
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
  error ('Matrices should have equal numbers of columns.');
end

ma  = mean(a');
a   = a - ma' * ones(1,ca);
mb  = mean(b');
b   = b - mb' * ones(1,cb);

aa = sum(a.*a,2);
bb = sum(b.*b,2)';
D  = (a*b') ./ sqrt(repmat(aa,1,rb) .* repmat(bb,ra,1));
D  = 0.5 * (1 - D);

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
