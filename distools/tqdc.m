%TQDC Trade-off Quadratic Discriminant (Regularized Bayes Normal Classifier)
%
%   W = TQDC(A,ALF,R,S,DIM)
%
% INPUT
%   A   NxK dataset (N points in a K-dimensional space)
%   ALF Trade-off parameter, ALF in [0,1] (optional; default: ALF = 0.1)
%   R,S Additional regularization parameters, 0 <= R,S <= 1
%       (optional; default: no regularization, i.e. R,S = 0)
%   DIM Dimension of subspace structure in covariance matrix (default: K)
%
% OUTPUT
%   W   Quadratic Bayes Normal Classifier mapping
%
% DESCRIPTION
% Computation of the quadratic classifier between the classes of the dataset
% A assuming normal densities. Each class covariance matrix Gi (i=1..C) is
% modeled as a convex combination between the original class covariance Gi and
% the diagonal marix Gdiag retrieved from the overall weighted (by priors)
% covariance matrix. So,
%     Gi = (1-ALF)*Gi + ALF*Gdiag
% If ALF=0, then you will get QDC.
% If ALF=1, then you will get NMSC.
%
% R and S (0 <= R,S <= 1) are additional parameters used for regularizing the
% resulting covariance matrices by
%     Gi = (1-R-S)*Gi + R*diag(diag(Gi)) + S*mean(diag(Gi))*eye(size(Gi,1))
% This covariance matrix is then decomposed as Gi = W*W' + sigma^2 * eye(K),
% where W is a KxM matrix containing the M leading principal components.
%
% The use of soft labels is supported. The classification A*W is computed by
% NORMAL_MAP.
%
% DEFAULT
% ALF = 0.1
% R   = 0
% S   = 0
% DIM = K (data dimension)
%
% EXAMPLES
% PREX_MCPLOT, PREX_PLOTC.
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd
% edition, John Wiley and Sons, New York, 2001.
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons,
% New York, 2002.
%
% SEE ALSO
% MAPPINGS, DATASETS, NMC, NMSC, LDC, UDC, QDC, QUADRC, NORMAL_MAP

% Copyright: R.P.W. Duin, E. Pekalska, D.M.J. Tax and P. Paclik
% ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function w = tqdc(a,alf,r,s,dim)

if (nargin < 5)
  prwarning(4,'Subspace dimensionality DIM not provided, assuming K.');
  dim = [];
end
if (nargin < 4)
  prwarning(4,'Regularisation parameter S not given, assuming 0.');
  s = 0;
end
if (nargin < 3)
  prwarning(4,'Regularisation parameter R not given, assuming 0.');
  r = 0;
end
if (nargin < 2)
  prwarning(4,'Trade-off parameter ALF not given, assuming 0.1.');
  alf = 0.1;
end

% No input arguments: return an UNTRAINED mapping
if (nargin < 1) | (isempty(a))
  w = mapping(mfilename,{alf,r,s,dim});
  w = setname(w,'Trade-off Bayes-Normal-2');
  return
end


% TRAIN the classifier
islabtype(a,'crisp','soft');
isvaldset(a,2,2);        % at least 2 objects per class and 2 classes

[m,k,c] = getsize(a);

% If the subspace dimensionality is not given, set it to the data dimensionality.
if (isempty(dim)),
  dim = k;
end

if (dim < 1) | (dim > k)
  error ('Number of dimensions DIM should lie in the range [1,K].');
end

% Assert whether A has the right labtype.
islabtype(a,'crisp','soft');

% Get mean vectors and class covariance matrices.
[U,G] = meancov(a);

% Calculate means and priors.
pars.mean  = +U;
pars.prior = getprior(a);

% in the NMSC limit:
Gtot = zeros(c,k);
for j = 1:c
  Gtot(j,:) = diag(G(:,:,j))';
end
Gtot = diag(pars.prior*Gtot);

% Calculate class covariance matrices.

pars.cov   = zeros(k,k,c);
for j = 1:c
  F = G(:,:,j);
  F = (1-alf)*F + alf*Gtot;

  % Regularize, if requested.
  if (s > 0) | (r > 0)
    F = (1-r-s) * F + r * diag(diag(F)) +s*mean(diag(F))*eye(size(F,1));
  end

  % If DIM < K, extract the first DIM principal components and estimate
  % the noise outside the subspace.

  if (dim < k)
    [eigvec,eigval] = eig(F);
    eigval = diag(eigval);
    [dummy,ind] = sort(-eigval);

    % Estimate sigma^2 as avg. eigenvalue outside subspace.
    sigma2 = mean(eigval(ind(dim+1:end)));

    % Subspace basis: first DIM eigenvectors * sqrt(eigenvalues).
    F = eigvec(:,ind(1:dim)) * diag(eigval(ind(1:dim))) * eigvec(:,ind(1:dim))' + sigma2 * eye(k);
  end
  pars.cov(:,:,j) = F;
end

w = mapping('normal_map','trained',pars,getlab(U),k,c);
w = setname(w,'Trade-off Bayes-Normal-2');
w = setcost(w,a);

return;
