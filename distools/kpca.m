%KPCA Kernel Principal Component Analysis
%
%   [W,L] = KPCA(K,ALF)
%         OR
%   [W,L] = KPCA(W,ALF)
%
% INPUT
%   K   NxN kernel or symmetric similarity matrix (dataset)
%   W   Trained KPCA projection a kernel-induced space
%   ALF Parameter determining the dimensionality and the mapping (optional, default: Inf)
%       (0,1) - fraction of the total preserved variance; e.g. ALF = 0.9
%       M     - Number of dimensions in the subspace; e.g. ALF = 5
%       Inf   - No dimensionality reduction, keeping all dimensions (it's VERY noisy)
%
% OUTPUT
%   W   Kernel PCA mapping
%   L   Sorted list of eigenvalues
%
% DEFAULTS
%   ALF = Inf
%
% DESCRIPTION
% Performs principal component analysis in a kernel-induced space by finding
% M significant directions. M is provided directly or determined by ALF such
% that the fraction ALF of the total variance is preserved. L is a sorted list
% of eigenvalues describing the variances in the kernel subspace.
% The projection X is found as X = K*W.
%
% A trained mapping can be reduced further by:
%   W = KPCA(W,ALF)
%
% SEE ALSO
% MAPPINGS, DATASETS, PCA, PSEM
%
% REFERENCE
% B. Scholkopf, A. Smola, and K.-R. Muller. Kernel Principal Component Analysis.
% in Advances in Kernel Methods - SV Learning, pages 327-352. MIT Press, Cambridge, MA, 1999.

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester, UK


function [W,L] = kpca (K,alf)

if nargin < 2 | isempty(alf),
  alf = inf;
end
if nargin < 1 | isempty(K),
  W = mapping(mfilename,alf);
  W = setname(W,'Kernel PCA');
  return
end


if (isdataset(K) | isa(K,'double')),
  if ismapping(alf),
    % APPLY MAPPING: project new data using the trained mapping.
    [m,n] = size(K);
    pars  = getdata(alf);

    % Parameters
    Q   = pars{1};  % Eigenvectors
    Kme = pars{2};  % Vector of averaged kernel values

    % Centering the kernel
    H = -repmat(1/n,n,n);
    H(1:n+1:end) = H(1:n+1:end) + 1;     % H = eye(n) - ones(n,n)/n
    K = (K - Kme(ones(m,1),:)) * H;
    W = K*Q;
    if isdataset(W),
      W.name = ['Projected '  updname(W.name)];
    end
    return;
  end
end


if ~isnumeric(alf) & ~isinf(alf)
  error('Wrong ALF.')
end
if alf <= 0
  error('ALF should be positive.')
end


% REDUCE ALREADY TRAINED MAPPING
if ismapping(K),
  pars  = getdata(K);

  Q  = pars{1};
  L  = pars{3};
  m  = size(Q,1);

  [ll,P] = sort(-abs(L));
  L = L(P);
  Q = Q(:,P);
  J = seleigs(L,alf);       % J is the index of selected eigenvalues
  Q = Q(:,J);               % Eigenvectors
  L = L(J);                 % Eigenvalues

  W = mapping(mfilename,'trained',{Q,pars{2},L},[],m,length(J));
  W = setname(W,'KPCA');
  return
end



% TRAIN MAPPING
K  = +K;
[n,m] = size(K);

% Tolerance value used in comparisons
if mean(diag(+K)) < 1,
  tol = 1e-12;
else
  tol = 1e-10;
end

if ~issym(K,tol),
  prwarning(1,'Kernel should be symmetric. It is made so by averaging.')
  K  = 0.5 * (K+K');
end

eigmin = min(eig(K));
if eigmin < 0,
  error(['K is not psd. Minimum eig(K) = ' ...
  num2str(eigmin) '. Please regularize the kernel appropriately or use IKPCA.']);
end

Kme = mean(K,2)';

% Project the data such that the mean lies at the origin
H = -repmat(1/n,n,n);
H(1:n+1:end) = H(1:n+1:end) + 1;  % H = eye(n) - ones(n,n)/n
K  = H * K * H;            % K is now the centered kernel
K  = 0.5 * (K+K');         % Make sure that K is symmetric after centering

[Q,L]  = eig(K);
Q      = real(Q);
l      = diag(real(L));
[lm,Z] = sort(-l);
Q      = Q(:,Z);
l      = l(Z);             % Eigenvalues are sorted by decreasing value


J = seleigs(l,alf);        % J is the index of selected eigenvalues
L = l(J);                  % Eigenvalues
Q = Q(:,J);                % Eigenvectors

% Normalize Q such that the eigenvectors of the covariance
% matrix are orthonormal
Q = Q* diag(1./sqrt(diag(Q'*K*Q)));

% Determine the mapping
W = mapping(mfilename,'trained',{Q,Kme,L},[],m,length(J));
W = setname(W,'KPCA');
return
