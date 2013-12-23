%BINDISTM Dissimilarity Matrix between Binary Vectors
%
%   D = BINDISTM(A,B,TYPE)
%
% INPUT
%   A     NxK Binary matrix or dataset
%   B     MxK Binary matrix or dataset (optional; default: B=A)
%   TYPE	Type of the similarity S (optional; default: 'J'):
%           'SM', 'Simple-match':   (a+d)/(a+b+c+d)
%           'RR', 'Russel-Rao':     a/(a+b+c+d)
%           'J',  'Jaccard':        a/(a+b+c)
%           'D',  'Dice':           a/(a+0.5*(b+c))
%           'SS', 'Sokal-Sneath':   (a+d)/(a+0.5*(b+c)+d)
%           'RT', 'Rogers-Tanimoto':(a+d)/(a+2*(b+c)+d)
%           'K',  'Kulczynski':     0.5*(a/(a+b) + a/(a+c))
%           'A1', 'Anderberg1':     a/(a+2*(b+c))
%           'A2', 'Anderberg2':     0.5*(a/(a+b) + a/(a+c) + d/(c+d) + d/(b+d))
%           'H',  'Hamman':         ((a+d)-(b+c))/(a+b+c+d)
%           'Y',  'Yule':           (a*d -b*c)/(a*d+b*c)
%           'P1', 'Pearson1':       (a*d)/sqrt((a+b)*(a+c)*(b+d)*(c+d))
%           'P2', 'Pearson2':       (a*d-b*c)/sqrt((a+b)*(a+c)*(b+d)*(c+d))
%           'O',  'Ochiai':         a/sqrt((a+b)*(a+c))
%         The distance D is computed as D=sqrt(1-S).
%
%         Type of distance:
%           'HG', 'Hamming':        (b+c)
%           'EU', 'Euclidean':      sqrt(b+c)
%           'VAR','Variance':       0.25*(b+c)/(a+b+c+d)
%           'BC', 'Bray-Curtis':    (b+c)/(2*a+b+c)
%           'SD', 'Size-diff':      (b-c)^2/(a+b+c+d)^2
%           'PD', 'Pattern-diff':   b*c/(a+b+c+d)^2
%           'SHD','Shape-diff':     ((a+b+c+d)*(b_c)-(b-c)^2)/(a+b+c+d)^2;
%
% OUTPUT
%   D 		NxM Dissimilarity matrix or dataset
%
% DESCRIPTION
% Distance between sets of binary vectors, A and B.
% The distances which are non-metric: 'K','A2','Y','BC','SD','PD' and 'SHD'.
% The distances which do not have a Euclidean behaviour: 'SS','K','A2','Y','HG',
% 'VAR','BC','SD','PD' and 'SHD'. If for the similarity S defined above, D is
% computed as D=1-S, then the following distances are non-metric: 'D','SS','K',
% 'A2','Y','P1','P2',and 'O', and all of the distances are non-Euclidean.
%
% NOTE
% In some cases the operations may be undefined such as 0/0. This results
% in NANs which are replaced here by zeros.
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% DEFAULT
%   B    = A
%   TYPE = 'J'
%
% REFERENCE
% J.Gower, Metric and Euclidean Properties od Dissimilarity Coefficients.
% Journal of Classification, no.5, 5-48, 1986.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function D = bindistm(A,B,type)

if nargin < 3,
  type = 'J';
end

bisa = (nargin < 2 | isempty(B));
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

if any(a~=0 & a~=1) | any(b~=0 & b~=1),
  error('Data should be binary.');
end

Aij = a*b';
Bij = a*(1-b)';
Cij = (1-a)*b';
Dij = (1-a)*(1-b)';

D = [];
switch lower(type)
  case {'hg','hamming'}
    D = (Bij+Cij);
  case {'eu','euclidean'}
    D = sqrt(Bij+Cij);
  case {'var','variance'}
    D = 0.25*(Bij+Cij)/ca;
  case {'bc','bray-curtis'}
    D = (Bij+Cij)./(2*Aij+Bij+Cij);
  case {'sd','size-diff'}
    D = (Bij-Cij).^2./ca^2;
  case {'pd','pattern-diff'}
    D = Bij.*Cij./ca^2;
  case {'shd','shape-diff'}
    D = (ca*(Bij_Cij)-(Bij-Cij).^2)./ca^2;
%
  case {'sm','simple-match'}
    S = (Aij+Dij) ./ ca;
  case {'rr','russel-rao'}
    S = Aij ./ ca;
  case {'j','jaccard'}
    S = Aij ./ (Aij+Bij+Cij);
  case {'d','dice'}
    S = Aij ./ (Aij+0.5*(Bij+Cij));
  case {'ss','sokal-sneath'}
    S = (Aij +Dij)./ (Aij + 0.5*(Bij+Cij) + Dij);
  case {'a1','anderberg1'}
    S = Aij./ (Aij + 2*(Bij+Cij));
  case {'rt','rogers-tanimoto'}
    S = (Aij +Dij)./ (Aij + 2*(Bij+Cij)+Dij);
  case {'k','kulczynski'}
    S = 0.5*(Aij./ (Aij + Bij) + Aij./ (Aij + Cij));
  case {'a2','anderberg2'}
    S = 0.5*(Aij./ (Aij + Bij) + Aij./ (Aij + Cij) + Dij./ (Cij + Dij) + Dij./ (Bij + Dij) );
  case {'h','hamman'}
    S = ((Aij + Dij) - (Bij + Cij))/ca;
  case {'y','yule'}
    S = (Aij .* Dij  - Bij .* Cij) ./ (Aij .* Dij  + Bij .* Cij);
  case {'p1','pearson1'}
    S = (Aij .* Dij) ./ sqrt((Aij + Bij) .* (Aij + Cij).*(Bij + Dij).*(Cij + Dij));
  case {'p2','pearson2'}
    S = (Aij .* Dij - Bij .* Cij) ./ sqrt((Aij + Bij) .* (Aij + Cij).*(Bij + Dij).*(Cij + Dij));
  case {'o','ochiai'}
    S = Aij / sqrt((Aij + Bij) .* (Aij + Cij));
  othwerwise
    error('Wrong type.');
end

if isempty(D),
  D = sqrt(1 - S);
end

% Replace potential NaNs by zeros
D(find(isnan(D))) = 0;

% Check numerical inaccuracy
D(find(D<eps)) = 0;

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
