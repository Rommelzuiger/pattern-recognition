%QDISTM  Distance Matrix for Quantitative Variables
%
%     D = QDISTM (A,B,TYPE,P)
%       OR
%     D = QDISTM (A,B)
%       OR
%     D = QDISTM (A,TYPE,P)
%       OR
%     D = QDISTM (A,TYPE)
%
% INPUT
%   A   NxK Matrix or dataset
%   B   MxK Matrix or dataset
%   TYPE  Type of the dissimilarity D (optional; default: 'E'):
%           'E',  'Euclidean'
%           'SQE','Square-Euclidean'
%           'LP', 'LP-distance'
%           'BC', 'Bray-Curtis'
%           'CAN','Canberra'
%           'COR','Correlation'
%           'COS','Cosine'
%           'DIV','Divergence'
%           'EXP','Exponent'
%           'S',  'Soergel'
%           'SAM','Spectral-Angular-Mapper'
%           'TAX','Taxonomic'
%           'WS', 'Ware-Hedges'
%   P   Parameter, P > 0 (optional, default: 1)
%
% OUTPUT
%   D   NxM Dissimilarity matrix or dataset
%
% DESCRIPTION
% Computation of the distance matrix D between two sets of vectors, A and B.
% Distances between vectors X and Y are computed as:
%   'E':   d(X,Y) = (sum_i (|X_i - Y_i|^2))^(1/2)
%   'SQE': d(X,Y) =  sum_i (|X_i - Y_i|^2)
%   'LP':  d(X,Y) = (sum_i (|X_i - Y_i|^P))^(1/P)
%   'BC':  d(X,Y) =  sum_i (|X_i - Y_i|)/sum_i (X_i + Y_i)
%   'CAN': d(X,Y) =  sum_i (|X_i - Y_i|)/sum_i (|X_i| + |Y_i|)
%   'COR': d(X,Y) =  (1 - COV(X,Y) / sqrt(Var(X) * VAR(Y)))/2
%   'COS': d(X,Y) =  (1 - X'*Y/(||X||*||Y||))
%   'DIV': d(X,Y) =  sum_i {|X_i - Y_i|^2/(X_i + Y_i)^2}
%   'EXP': d(X,Y) =  1 - exp (-(X-Y)'(X-Y)/P^2)%
%   'S':   d(X,Y) =  sum_i (|X_i - Y_i|)/max_i {X_i,Y_i}
%   'SAM': d(X,Y) =  P arcos (X'Y/P^2)
%   'TAX': d(X,Y) = (sum_i |X_i - Y_i|^P/r_i^P)^(1/P)
%   'WS':  d(X,Y) =  sum_i {1 - min_i{X_i,Y_i}/max_i{X_i,Y_i}}
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% DEFAULT
%   B    = A
%   TYPE = 'E'
%   P    = 1
%
% SEE ALSO
%   CORRDISTM, COSDISTM, DISTM, EXPDISTM, EUDISTM, LPDISTM, SAMDISTM,
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function D = qdistm (A,B,type,p)

bisa = 0;
if nargin < 2,
  p    = 1;
  B    = A;
  type = 'E';
  bisa = 1;
elseif nargin < 3,
  if isstr(B),
    p    = 1;
    type = B;
    B    = A;
    bisa = 1;
  else
    p    = 1;
    type = 'E';
  end
elseif nargin < 4,
  if ~isstr(type),
    p    = type;
    type = B;
    B    = A;
    bisa = 1;
  else
    p    = 1;
  end
else
  ;
end

if ~isstr(type)
  error ('TYPE is a string.');
end

if p <= 0,
  error ('The parameter P must be positive.');
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
switch lower(type)
  case {'e','euclidean'}
    D = sqrt(distm(a,b));
  case {'sqe','square-euclidean'}
    D = distm(a,b);
  case {'lp','lp-distance'}
    D = lpdistm(a,b,p);
  case {'bc','bray-curtis'}
    for i=1:rb
      D(:,i) = sum(abs(repmat(b(i,:),ra,1) - a),2);
      D(:,i) = D(:,i) ./ sum((repmat(b(i,:),ra,1) + a),2);
    end
  case {'can','canberra'}
    for i=1:rb
      D(:,i) = sum( abs(repmat(b(i,:),ra,1) - a) ./ (repmat(abs(b(i,:)),ra,1) + abs(a)), 2);
    end
  case {'cor','correlation'}
    D = corrdistm(a,b);
  case {'cos','cosine'}
    D = cosdistm(a,b);
  case {'div','divergence'}
    for i=1:rb
      Z = (abs(repmat(b(i,:),ra,1) - a)).^p;
      D(:,i) = sum (Z ./(repmat(b(i,:),ra,1) + a).^p, 2);
      D(:,i) = D(:,i).^(1/p);
      clear Z;
    end
  case {'exp','exponent'}
    D = expdistm(a,b,p);
  case {'sam','spectral-angular-mapper'}
    D = samdistm(a,b,p);
  case {'s','soergel'}
    for i=1:rb
      D(:,i) = sum(abs(repmat(b(i,:),ra,1) - a),2) ./ sum(max(repmat(b(i,:),ra,1),a),2);
    end
  case {'tax','taxonomic'}
    rr = max(b) - min(b);
    for i=1:rb
      D(:,i) = sum( (abs(repmat(b(i,:),ra,1) - a)./repmat(rr,ra,1)).^p,2);
      D(:,i) = D(:,i).^(1/p);
    end
  case {'ws','ware-hedges'}
    for i=1:rb
      D(:,i) = sum(1 -  min(repmat(b(i,:),ra,1),a) ./ max(repmat(b(i,:),ra,1),a),2);
    end
  otherwise
    error('Wrong dissimilarity type.');
end

if bisa,
  D = 0.5*(D+D');         % Make sure that distances are symmetric for D(A,A)
end

% Set object labels and feature labels
if xor(isda, isdb),
  prwarning(1,'One matrix is a dataset and the other is not. ')
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
