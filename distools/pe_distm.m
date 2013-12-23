%PE_DISTM Square Pseudo-Euclidean (PE) Distance Between Two Datasets
%
%   D = PE_DISTM(A)
%     OR
%   D = PE_DISTM(A,B)
%
% INPUT
%   A    PE dataset of size NxK
%   B    PE dataset of size MxK (default A = B)
%
% OUTPUT
%   D   NxM dissimilarity matrix or dataset
%
% DESCRIPTION
% Computation of the square pseudo-Euclidean distance matrix D between two sets
% of vectors, A and B. The pseudo-Euclidean distance with the signature SIG
% (e.g. SIG = [10 5]) between vectors X and Y is computed as an indefinite
% 'Euclidean' distance:
%     D(X,Y) = (X-Y)'*J*(X-Y),
% where J is a diagonal matrix with 1's, followed by -1's.
% J = diag ([ONES(SIG(1),1);  -ONES(SIG(2),1)]);
%
% In a PE dataset the signature is stored in the user field, see
% SETSIG. This signature is derived from A. It is not stored in D as D
% does not contain vectors in a PE space.
%
% D is a dataset with the labels defined by the labels of A and feature labels 
% defined by the labels of B. 
%
% REMARKS
% Note that square pseudo-Euclidean distances can be negative.
%
% SEE ALSO
% DATASET, SETSIG, DISTM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function D = pe_distm(A,B)

  isdataset(A);
  bisa = 0;
  if nargin < 2, B = A; bisa = 1; end
  sig = getsig(A);

  if ~isdataset(B), B = dataset(B,1); end

  a    = +A;
  b    = +B;
  [ra,ca] = size(a);
  [rb,cb] = size(b);

  if ca ~= cb,
    error ('The datasets should have the same number of features.');
  end

  if any(sig) < 0 | sum(sig) ~= ca,
    error('Signature vector SIG is invalid: its sum should equal feature size');
  end

  J = [ones(1,sig(1))  -ones(1,sig(2))];
  D = - 2 .* a * diag(J) * b';
  D = D + ones(ra,1) * (J*(b'.*b'));
  D = D + (J * (a'.*a'))' * ones(1,rb);

  if bisa
    D = 0.5*(D+D');         % Make sure that distances are symmetric for D(A,A)
  end

  % Set object labels and feature labels
  D = setdata(A,D,getlab(B));
  D.name = 'Square Pseudo-Euclidean distance matrix';
  if ~isempty(A.name)
    D.name = [D.name ' for ' A.name];
  end

return
