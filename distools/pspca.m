%PSPCA Pseudo-Euclidean Principal Component Analysis
%
%     [W,SIG,L] = PSPCA(X,XSIG,ALF)
%
% INPUT
%   X     NxK data
%   XSIG  Signature of the input pseudo-Euclidean space; (default: [K 0])
%   ALF   Parameter determining the dimensionality and the mapping (optional, default: Inf)
%         (0,1)   - fraction of the total (absolute value) preserved variance
%         Inf     - no dimensionality reduction, keeping all dimensions (it's noisy)
%         'p'     - projection into a Euclidean space based on positive eigenvalues only
%         'PARp'  - projection into a Euclidean space based on the PAR fraction of
%                   positive eigenvalues; e.g. ALF = '0.9p'
%         'n'     - projection into a Euclidean space based on negative eigenvalues only
%         'PARn'  - projection into a (negative) Euclidean space based on the PAR fraction
%                   of negative eigenvalues; e.g. ALF = '0.7n'
%         'P1pP2n'- projection into a Euclidean space based on the P1 positive eigenvalues
%                   and P2 negative eigenvalues; e.g. ALF = '0.7p0.1n', ALF = '7p2n'
%         1 .. N  - number of dimensions in total
%         [P1 P2] - P1 dimensions or preserved fraction of variance in the positive subspace
%                   and P2 dimensions or preserved fraction of variance in the negative
%                   subspace; e.g. ALF = [5 10], ALF = [0.9 0.1]
%
% OUTPUT
%   W     PCA mapping in a pseudo-Euclidean space
%   SIG   Signature of the output pseudo-Euclidean space
%   L     List of eigenvalues
%
% DEFAULT
%   XSIG = [K 0]
%   ALF  = INF
%
% DESCRIPTION
% PCA mapping W from a K-dimensional pseudo-Euclidean space to an M-dimensional
% pseudo-Euclidean subspace. M is determined by ALF. The subspace is found, e.g.
% such that at least a fraction ALF of the total variance is preserved for ALF
% in (0,1). The resulting Y is found by X*W. The parameter SIG describes the
% signature of the subspace. L is a sorted list of eigenvalues describing the
% variances in the (pseudo-)Euclidean space.
%
% If X is a labeled dataset, then the averaged covariance matrix is weighted
% by class priors.
%
% Note that a PCA decomposition in a pseudo-Euclidean space is different than in
% a Euclidean space. Namely, CJ = Q*L*inv(Q), where CJ is a pseudo-Euclidean
% covariance matrix computed such that CJ= C*J, where C is a Euclidean covariance
% matrix, J is the fundamental symmetry (taking part in inner products). Q is
% J-orthogonal, i.e. Q'*J*Q = J, hence inv(Q) = J*Q'*J.
%
% SEE ALSO
% MAPPINGS, DATASETS, PCA, KPSEM, PSEM
%
% LITERATURE
% 1. E. Pekalska, R.P.W. Duin, The Dissimilarity representation in Pattern Recognition.
%    Foundations and Applications. World Scientific, Singapore, 2005.
% 2. L. Goldfarb, A unified approach to pattern recognition, Pattern Recognition, vol.17, 575-582, 1984.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function [W,outsig,L,Q] = pspca(a,sig,alf,prec)

if nargin < 4 | isempty(prec), prec = 1e-4; end
if nargin < 3 | isempty(alf), alf = inf; end
if nargin < 2 | isempty(sig), sig = [size(a,1) 0]; end
if nargin < 1 | isempty(a),
  W = mapping(mfilename,sig,alf,prec);
  W = setname(W,'Pseudo-Euclidean PCA');
  return
end


if (isdataset(a) | isa(a,'double')),

  if ismapping(sig),
    % APPLY MAPPING: project new data using the trained mapping.
    [m,n] = size(a);
    pars  = getdata(sig);

    % Parameters
    v   = pars{1};  % Mapping that shifts data to the origin
    JQ  = pars{2};  % J*Q
    sig = pars{3};  % Signature in the output space
    W = (a*v) * JQ;
    if isdataset(W),
      W.user = sig;
      W.name = updname(W.name);
    end
    return;
  end
end


% TRAIN THE MAPPING
[m,k] = size(a);
if m < 2,
  error('At least two objects are expected.');
end
if sum(sig) ~= k,
  error('Signature does not fit the data dimensionality.')
end
isdset = isdataset(a);


% Shift mean of data to the origin
v  = scalem(+a);
aa = a*v;

if ~isdset,    % Unlabeled data
  A = +aa;
else
  c = max(getnlab(aa));
  if c == 0,
    A = +aa;
  else
    p = getprior(a);
    A = [];
    for j = 1:c
      A = [A; +seldat(aa,j)*p(j)];
    end
  end
end
G = cov(A);
G = 0.5*(G+G');     % Make sure G is symmetric
if sig(2) > 0,
  J = diag([ones(sig(1),1); -ones(sig(2),1)]);
  G = G*J;
end

[Q,L]  = eig(G);
Q      = real(Q);
l      = diag(real(L));
[lm,Z] = sort(-abs(l));
Q      = Q(:,Z);
l      = l(Z);

[I,outsig] = seleigs(l,alf,prec); % I is the index of selected eigenvalues
L = l(I);                         % Eigenvalues
Q = Q(:,I);                       % Eigenvectors


if sig(2) > 0,
  % Q is NOT orthogonal, but should be J-orthogonal, i.e. Q'*J*Q = J
  % Normalize Q to be J-orthogonal
  Q = Q*diag(1./sqrt(abs(diag(Q'*J*Q))));
  Q = J*Q;
end


% Determine the mapping
W = mapping(mfilename,'trained',{v,Q,outsig,sig},[],k,sum(outsig));
W = setname(W,'Pseudo-Euclidean PCA');
return
