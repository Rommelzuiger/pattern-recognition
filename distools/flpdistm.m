%FLPDISTM  lp (p > 0) (Non)-Metric Distance Matrix
%
%     D = FLPDISTM (A,B,P)
%       OR
%     D = FLPDISTM (A,B)
%       OR
%     D = FLPDISTM (A,P)
%       OR
%     D = FLPDISTM (A)
%
% INPUT
%   A   NxK Matrix or dataset
%   B   MxK Matrix or dataset
%   P   Parameter; P > 0
%
% OUTPUT
%   D   NxM Dissimilarity matrix or dataset
%
% DEFAULT
%   P = 1
%   B = A
%
% DESCRIPTION
% Fast computation of the distance matrix D between two sets of vectors, A and B.
% This can ONLY be used for small sets A and B as the memory is significantly used
% by computing 3D matrices of the size M x N x K.
% Distances between vectors X and Y are computed using the lp distance:
%     d(X,Y) = (sum (|X_i - Y_i|.^P))^(1/P)
%                i
% If P = Inf, then the max norm distance is computed:
%     d(X,Y) = max (|X_i - Y_i|)
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% DEFAULT
%   P = 1
%   B = A
%
% REMARKS
%   P >= 1      => D is metric
%   P in (0,1)  => D is non-metric; D.^P is metric and l1-embeddable
%   P = 1/2     => D is city block / Euclidean distance
%
% SEE ALSO
%   LPDISTM, EUDISTM
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function D = flpdistm (A,B,p)
bisa = 0;
if nargin < 2,
  p = 1;
  B = A;
  bisa = 1;
else
  if nargin < 3,
    if max (size(B)) == 1,
      p = B;
      bisa = 1;
      B = A;
    else
      p = 1;
    end
  end
end

if p <= 0,
  error ('The parameter p must be positive.');
end

isda = isdataset(A);
isdb = isdataset(B);
a    = +A;
b    = +B;
[ra,ca] = size(a);
[rb,cb] = size(b);

if ca ~= cb,
  error ('The matrices should have the same number of columns.');
end

D = zeros(ra,rb);
if p < Inf,
  D = sum ((abs (repmat(permute(a,[1 3 2]), [1 rb 1]) - ...
                 repmat(permute(b,[3 1 2]), [ra 1 1]))).^p,3).^(1/p);
else
  D = max ((abs (repmat(permute(a,[1 3 2]), [1 rb 1]) - ...
                 repmat(permute(b,[3 1 2]), [ra 1 1]))),[],3);
end


% Check numerical inaccuracy
D (find (D < eps)) = 0;   % Make sure that distances are nonnegative
if bisa,
  D = 0.5*(D+D');         % Make sure that distances are symmetric for D(A,A)
end


if xor(isda, isdb),
  prwarning(1,'One matrix is a dataset and the other not. The result is a matrix.')
elseif isda & isdb,
  D = setdata(A,D,getlab(B));
  D.name = 'Distance matrix';
  if ~isempty(A.name)
    D.name = [D.name ' for ' A.name];
  end
else
  ;
end
return
