%SAMDISTM  Distance matrix based on Spectral Angular Mapper (SAM)
%          distance, which is also the spherical geodesic distance
%
%   D = SAMDISTM (A,B,R)
%     OR
%   D = SAMDISTM (A,B)
%     OR
%   D = SAMDISTM (A,R)
%     OR
%   D = SAMDISTM (A)
%
% INPUT
%   A   NxK matrix (dataset)
%   B   MxK matrix (dataset)
%   R   Radius (optional, default: 1)
%
% OUTPUT
%   D   NxM dissimilarity matrix (dataset)
%
% DESCRIPTION
% Computes the distance matrix D between two sets of vectors, A and B.
% Distances between vectors X and Y are computed based on the spherical
% geodesic formula:
%     D(X,Y) = R arcos (X'Y/R^2)
% X and Y are normalized to a unit length.
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% DEFAULT
%   R = 1
%
% REMARKS
% A square SAM-distance D(A,A) for a finite set A can be proved to be
% the l_1-distance (LPDISTM). D(A,A).^{1/2} has a Euclidean behavior, so
% it can be embedded by PSEM in a Euclidean space.
%

% SEE ALSO
% JACSIMDISTM, CORRDISTM, LPDISTM, DISTM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function D = samdistm (A,B,r)

bisa = 0;
if nargin < 2,
  r = 1;
  B = A;
  bisa = 1;
else
  if nargin < 3,
    if max (size(B)) == 1,
      r = B;
      B = A;
      bisa = 1;
    else
      r = 1;
    end
  end
end

if r <= 0,
  error ('The parameter R must be positive.');
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

aa = sum(a.*a,2);
bb = sum(b.*b,2)';
D  = (a*b') ./sqrt(aa(:,ones(rb,1)) .* bb(ones(ra,1),:));
D  = r * acos(D/r^2);

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
