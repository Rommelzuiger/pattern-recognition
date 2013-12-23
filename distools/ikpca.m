%IKPCA Indefinite Kernel Principal Component Analysis
%
%   [W,SIG,L] = IKPCA(K,ALF)
%             OR
%   [W,SIG,L] = IKPCA(W,ALF)
%
% INPUT
%   K   NxN kernel or symmetric similarity matrix (dataset)
%   W   Trained IKPCA projection into a kernel-induced pseudo-Euclidean subspace
%   ALF Parameter determining the dimension and the mapping (optional, default: Inf)
%       (0,1) - fraction of the total magnitude of the preserved variance
%       Inf   - no dimension reduction; keeping all dimensions (it's VERY noisy)
%       'p'   - projection into a Euclidean space based on the positive eigenvalues only
%       'PARp'- projection into a Euclidean space based on the PAR fraction of
%               positive eigenvalues; e.g. ALF = '0.9p'
%       'n'   - projection into a Euclidean space based on the negative eigenvalues only
%       'PARn'- projection into a (negative) Euclidean space based on the PAR fraction
%               of negative eigenvalues; e.g. ALF = '0.7n'
%     'P1pP2n'- projection into a Euclidean space based on the P1 positive eigenvalues
%               and P2 negative eigenvalues; e.g. ALF = '0.7p0.1n', ALF = '7p2n'
%     1 .. N  - total number of significant dimensions
%     [P1 P2] - P1 dimensions or a fraction of preserved variance in the positive
%               subspace and P2 dimensions or a fraction of preserved variance in
%               the negative subspace; e.g. ALF = [5 10], ALF = [0.9 0.1]
%
% OUTPUT
%   W   Indefinite Kernel PCA mapping
%   SIG Signature of the kernel Pseudo-Euclidean subspace
%   L   Sorted list of eigenvalues
%
% DEFAULTS
%   ALF = Inf
%
% DESCRIPTION
% Performs principal component analysis in a kernel space (either Hilbert or Krein
% space) by finding M significant principal components. M is determined by ALF,
% e.g. such that at least the fraction ALF of the total variance is preserved.
% The parameter SIG describes the signature of the subspace. L is a sorted list
% of eigenvalues describing the variances in the (pseudo-)Euclidean subspace.
% The projection X is found as X = K*W. The signature is also stored in X.user.
%
% A trained mapping can be reduced by:
%   [W,SIG,L] = IKPCA(W,ALF)
%
% SEE ALSO
% MAPPINGS, DATASETS, PCA, PSEM
%
% REFERENCE
% 1. E. Pekalska and R.P.W. Duin. Indefinite Kernel Principal Component Analysis.
%    Technical Report, 2006.
% 2. B. Schölkopf, A. Smola, and K.-R. Müller. Kernel Principal Component Analysis.
%    in Advances in Kernel Methods - SV Learning, pages 327-352. MIT Press, Cambridge, MA, 1999.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester, UK


function [W,sig,L] = ikpca (K,alf,prec)

% PREC is the precision parameter used for the automatic
% selection (heuristic) of the number of dominant eigenvalues.
% This happens when SELEIGS is called with the parameter 'CUT'.

if nargin < 3 | isempty(prec),
  prec = 1e-4;
end
if nargin < 2 | isempty(alf),
  alf = inf;
end
if nargin < 1 | isempty(K),
  W = mapping(mfilename,alf,prec);
  W = setname(W,'IKPCA');
  return
end


if (isdataset(K) | isa(K,'double')),
  if ismapping(alf),
    % APPLY THE MAPPING: project new data using the trained mapping
    [m,n] = size(K);
    pars  = getdata(alf);

    % Parameters
    QL  = pars{1};    % Normalized eigenvectors
    Kme = pars{2};    % Row vector of the average original kernel values

    % Centering the kernel
    H = -repmat(1/n,n,n);
    H(1:n+1:end) = H(1:n+1:end) + 1;     % H = eye(n) - ones(n,n)/n
    K = (K - repmat(Kme,m,1)) * H;
    W = K*QL;
    if isdataset(W),
      W.name = ['Projected '  updname(W.name)];
      W.user = pars{4};  % Signature of the PE subspace
    end
    return;
  end
end



% REDUCE ALREADY TRAINED MAPPING
if ismapping(K),
  pars  = getdata(K);

  QL = pars{1};
  L  = pars{3};
  m  = size(Q,1);

  [ll,P]  = sort(-abs(L));
  L       = L(P);
  QL      = QL(:,P);
  [J,sig] = seleigs(L,alf,pars{5});% J is the index of selected eigenvalues
  QL      = QL(:,J);               % Eigenvectors
  L       = L(J);                  % Eigenvalues

  W = mapping(mfilename,'trained',{QL,pars{2},L,sig,pars{5}},[],m,length(J));
  W = setname(W,'IKPCA');
  return
end




% TRAIN THE MAPPING
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
if eigmin < 0 & abs(eigmin) < 1e-12
  disp(['Minimum eig(K) = ' num2str(eigmin) '. K might be indefinite due to numerical inaccuracies.']);
end


Kme = mean(K,2)';

% Project the data such that the mean lies at the origin
H = -repmat(1/n,n,n);
H(1:n+1:end) = H(1:n+1:end) + 1;  % H = eye(n) - ones(n,n)/n
K  = H * K * H;            % K is now a centered kernel
K  = 0.5 * (K+K');         % Make sure that K is symmetric after centering

[Q, L] = eig(K);
Q      = real(Q);
l      = diag(real(L));
[lm,Z] = sort(-abs(l));
Q      = Q(:,Z);
l      = l(Z);

% Eigenvalues are now sorted by decreasing absolute value

[J,sig] = seleigs(l,alf,prec);  % J is the index of selected eigenvalues
L = l(J);                       % Eigenvalues
Q = Q(:,J);                     % Eigenvectors

% Normalize Q such that the eigenvectors of the corresponding
% covariance matrix are J-orthonormal
QL = Q * diag(1./sqrt(abs(L)));

% Determine the mapping
W = mapping(mfilename,'trained',{QL,Kme,L,sig,prec},[],m,sum(sig));
W = setname(W,'IKPCA');
return
