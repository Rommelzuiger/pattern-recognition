%PSDISTM Square Pseudo-Euclidean Distance Between Two Datasets
%
%   D = PSDISTM(A,SIG)
%     OR
%   D = PSDISTM(A,B,SIG)
%
% INPUT
%   A    NxK Matrix or dataset
%   B    MxK Matrix or dataset
%   SIG  2x1 vector being the signature of the pseudo-Euclidean space;
%        SIG(1)+SIG(2)=K
%
% OUTPUT
%   D   NxM dissimilarity matrix or dataset
%
% DESCRIPTION
% Computation of the square pseudo-Euclidean distance matrix D between two sets
% of vectors, A and B. The pseudo-Euclidean distance with the signature SIG
% (e.g. SIG = [10 5]) between vectors X and Y is comuted as an indefinite
% 'Euclidean' distance:
%     D(X,Y) = (X-Y)'*J*(X-Y),
% where J is a diagonal matrix with 1, followed by -1.
% J = diag ([ONES(SIG(1),1);  -ONES(sig(2),1)]);
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.
%
% REMARKS
% Note that square pseudo-Euclidean distances can be negative.
%
% SEE ALSO
% DISTM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function D = psdistm(A,B,sig)

bisa = 0;
if nargin < 2,
  error ('Inputs not specified');
elseif nargin < 3
  if max (size(B)) == 2 & min(size(B)) == 1,
    sig = B;
    B = A;
    bisa = 1;
  else
    error('Signature vector SIG expected.');
  end
else
  ;
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

if any(sig) < 0 | sum(sig) ~= ca,
  error('Signature vector SIG is invalid.');
end

J = [ones(1,sig(1))  -ones(1,sig(2))];
D = - 2 .* a * diag(J) * b';
D = D + ones(ra,1) * (J*(b'.*b'));
D = D + (J * (a'.*a'))' * ones(1,rb);

% Check numerical inaccuracy
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
  D.name = 'Square Pseudo-Euclidean distance matrix';
  if ~isempty(A.name)
    D.name = [D.name ' for ' A.name];
  end
end
return
