%RANKDISTM  Distance matrix between two data sets based on ranking
%
%     D = RANKDISTM (A,B,P)
%       or
%     D = RANKDISTM (A,B)
%       or
%     D = RANKDISTM (A,P)
%       or
%     D = RANKDISTM (A)
%
% INPUT
%   A   NxK Matrix or dataset
%   B   MxK Matrix or dataset
%   P   Parameter:
%       Integer - 1 .. K or 'MIN', 'MAX', 'MEDIAN' (optional; default: 'MEDIAN')
%
% OUTPUT
%   D   NxM dissimilarity matrix or dataset
%
% DESCRIPTION
% Computes the distance matrix D between two sets of vectors, A and B.
% Given the vectors X and Y, distances are computed using the ranked
% distance as:
%   D(X,Y) = P-th value of (sort {|X_1 - Y_1|, |X_2 - Y_2|,..,|X_K - Y_K|})
%
% For instance, for P = 1, the ranked distance becomes the minimum value of
% the differences |X_i - Y_i|, or for P = K, the infinty norm.
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% DEFAULT
% P = 'MEDIAN'
%
% SEE ALSO
% LPDISTM, EUDISTM, SIMDISTM, JACSIMDISTM, CORRDISTM, COSDISTM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function D = rankdistm (A,B,kk)
bisa = 0;
if nargin == 3,
  k = whichk(kk,ca);
elseif nargin < 2,
  k = 0;                  % median
  B = A;
  bisa = 1;
else
  k = whichk(B,ca);
  B = A;
  bisa = 1;
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
switch k,
  case 0,
    for i=1:rb
      D(:,i) = median (abs(repmat(b(i,:),ra,1) - a),2);
    end
    % D = median((abs (repmat (permute(a,[1 3 2]), [1 rb 1]) - ...
    %                  repmat (permute(b,[3 1 2]), [ra 1 1]))),3);
  case 1,
    for i=1:rb
      D(:,i) = min (abs(repmat(b(i,:),ra,1) - a),[],2);
    end
%    D = min((abs (repmat (permute(a,[1 3 2]), [1 rb 1]) - ...
%                  repmat (permute(b,[3 1 2]), [ra 1 1]))),[],3);
  case ra,
    for i=1:rb
      D(:,i) = max (abs(repmat(b(i,:),ra,1) - a),[],2);
    end
%    D = max((abs (repmat (permute(a,[1 3 2]), [1 rb 1]) - ...
%                  repmat (permute(b,[3 1 2]), [ra 1 1]))),[],3);
  otherwise
    for i=1:rb
      aa = sort (abs(repmat(b(i,:),ra,1) - a),2);
      D(:,i) = aa(:,k);
    end
%    aa = sort (abs (repmat (permute(a,[1 3 2]), [1 rb 1]) - ...
%               repmat (permute(b,[3 1 2]), [ra 1 1])), 3);
%    D  = aa(:,:,k);
end

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




function k = whichk(kk,ca)
if isstr(kk),
  switch lower(kk),
    case 'min'
      k = 1;
    case 'max'
      k = ca;
    case 'median',
      k = 0;
    otherwise
      error ('Wrong parameter k.');
  end
elseif max(size(kk)) == 1,
  k = kk;
else
  error ('Wrong parameter k.');
end

if k < 0 | k > ca,
  error ('The parameter k, if an integer, must be positive and not larger then the number of features.');
end
