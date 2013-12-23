%PE_EM Pseudo-Euclidean Linear Embedding
%
%     W = PE_EM(D,ALF,P)
%                 OR
%     W = PE_EM(W,ALF)
%
% INPUT
%   D   NxN symmetric dissimilarity matrix (dataset)
%   W   Trained linear embedding into a pseudo-Euclidean space
%   ALF Parameter determining the dimensionality and the mapping (optional, default: Inf)
%       (0,1)   - Fraction of the total (absolute value) preserved variance
%        Inf    - No dimensionality reduction, keeping all dimensions (it's VERY noisy)
%       'p'     - Projection into a Euclidean space based on positive eigenvalues only
%       'PARp'  - Projection into a Euclidean space based on the PAR fraction of
%                 positive eigenvalues; e.g. ALF = '0.9p'
%       'n'     - Projection into a Euclidean space based on negative eigenvalues only
%       'PARn'  - Projection into a (negative) Euclidean space based on the PAR fraction
%                 of negative eigenvalues; e.g. ALF = '0.7n'
%       'P1pP2n'- Projection into a Euclidean space based on the P1 positive eigenvalues
%                 and P2 negative eigenvalues; e.g. ALF = '0.7p0.1n', ALF = '7p2n'
%       1 .. N  - Number of dimensions in total
%       [P1 P2] - P1 dimensions or preserved fraction of variance in the positive subspace
%                 and P2 dimensions or preserved fraction of variance in the negative
%                 subspace; e.g. ALF = [5 10], ALF = [0.9 0.1]
%   P   Integer between 0 and N specifying which object is mapped at the origin;
%       0 stands for the mean; (optional, default: 0)
%
% OUTPUT
%   W   Linear embedding into a pseudo-Euclidean space
%
% DEFAULT
%   P   = 0
%   ALF = INF
%
% DESCRIPTION
% Linear mapping W onto an M-dimensional Pseudo-Euclidean _PE) subspace from a 
% symmetric, square dissimilarity matrix D such that the dissimilarities are preserved.
% M is determined by ALF. E.g., the subspace is found such that at least a fraction
% ALF of the total variance is preserved for ALF in (0,1). The resulting X is found
% by D*W. The signature of the obtained PE space (numbers of positive and negative
% directions) can be found by PE_SIG(W). The spectrum of the obtained space
% can be found by PE_SPEC(W). 
%
% A trained mapping can be reduced further by:   W = PE_EM(W,ALF)
% The signature of the obtained PE space can be found by PE_SIG(W)
% The spectrum of 
%
% SEE ALSO
% MAPPINGS, DATASETS, AUGPSEM, PCA, PE_PCA, PE_SPEC, GETSIG, SETSIG
%
% LITERATURE
% 1. L. Goldfarb, A unified approach to pattern recognition, Pattern Recognition, vol.17,
% 575-582, 1984.
% 2. E. Pekalska, P. Paclik, and R.P.W. Duin, A Generalized Kernel Approach to
% Dissimilarity-based Classification, Journal of Machine Learning Research,
% vol.2, no.2, 175-211, 2002.

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester

% This routine is adapted from PSEM. 
% Old PSEM mappings W may be converted by W = PE_EM(W)

function [W,sig,L,Q] = pe_em(d,alf,pzero,prec)

% PREC is the precision parameter used for the automatic
% selection (heuristic) of the number of dominant eigenvalues.
% This happens when SELEIGS is called with the parameter 'CUT'.

if nargin < 4, prec = [];  end
if nargin < 3, pzero = []; end
if nargin < 2  alf = [];   end

if nargin < 1 | isempty(d),
  W = mapping(mfilename,{alf,pzero,prec});
  W = setname(W,'PE embedding');
  return
end

if isempty(prec),  prec = 1e-4; end
if isempty(pzero), pzero = 0;   end
if isempty(alf),   alf = inf;   end

if (isdataset(d) | isa(d,'double'))
  if ismapping(alf)
    % APPLY THE MAPPING
    [m,n] = size(d);
    d     = d.^2;

    Q   = getdata(alf,'evec');  % Eigenvectors
    me  = getdata(alf,'mdis');  % Vector of the average squared original dissimilarities
    p   = getdata(alf,'mean');  % p=0 -> the mean of the embedded configuration lies at 0, 
                                % otherwise, it lies at pzero
    L   = getdata(alf,'eval');  % Eigenvalues

    % Project new data depending on p
    % (whether the mean or other object lies at the origin)
    if p == 0,
      H = -repmat(1,n,n)/n;
      H(1:n+1:end) = H(1:n+1:end) + 1;     % H = eye(n) - ones(n,n)/n
      W = -0.5 * (d - me(ones(m,1),:)) * H * Q * diag(sqrt(abs(L))./L);
    else
      W =  0.5 * (d(:,p) * ones(1,n) + me(ones(m,1),:) - d) * Q * diag(sqrt(abs(L))./L);
    end

    % Store signature in the USER field
    if isdataset(W),
      W = setname(W,['Projected ' updname(W.name)]);
      W = setsig(W,getdata(alf,'sig'));
    end
    return
  end
end



% REDUCE A TRAINED MAPPING
if ismapping(d)
  data = getdata(d);
  if iscell(data) % check for old type of mapping
    ispsem(d);
    dat.evec = data{1};
    dat.mdis = data{2};
    dat.mean = data{3};
    dat.eval = data{4};
    dat.sig  = data{5};
    dat.prec = data{6};
    d = setmapping_file(d,mfilename);
    d = setdata(d,dat);
    if nargin < 2
      W = d;
      return % conversion only
    else 
      data = dat;
    end
  end
  Q   = data.evec;  % Eigenvectors
  L   = data.eval;  % Eigenvalues
  m    = size(Q,1);

  [ll,K] = sort(-abs(L));
  L = L(K);
  Q = Q(:,K);
  [J,sig] = seleigs(L,alf,getdata(d,'prec'));
  data.evec = Q(:,J);        % Eigenvectors
  data.eval = L(J);          % Eigenvalues
 
  W = mapping(mfilename,'trained',data,[],m,length(J));
  W = setname(W,'PE embedding');
  return
end



% TRAIN THE MAPPING
% Tolerance value used in comparisons
if mean(+d(:)) < 1,
  tol = 1e-12;
else
  tol = 1e-10;
end

[n,m] = size(d);
if ~issym(d,tol),
  prwarning(1,'Matrix should be symmetric. It is made symmetric by averaging.')
  d = 0.5*(d+d');
end

if pzero > n,
  error('Wrong third parameter.');
end

d = (+d).^2;

if pzero == 0,
  % Project the data such that the mean lies at the origin
  H = -repmat(1/n,n,n);
  H(1:n+1:end) = H(1:n+1:end) + 1; % H = eye(n) - ones(n,n)/n
  H = -0.5 * H * d * H;            % H is now the matrix of inner products
else
  % Project the data such that pzero's object lies at the origin
  H = 0.5 * (d(:,pzero) * ones(1,n) + ones(n,1) * d(:,pzero)' - d);
end
H = 0.5*(H+H');                    % Make sure H is symmetric

[Q,L]  = eig(H);
Q      = real(Q);
l      = diag(real(L));
[lm,Z] = sort(-abs(l));
Q      = Q(:,Z);
l      = l(Z);                     % Eigenvalues are sorted by decreasing absolute value

[J,sig] = seleigs(l,alf,prec);     % J is the index of the selected eigenvalues
data.eval = l(J);                          % Eigenvalues
data.evec = Q(:,J);                        % Eigenvectors
data.mdis = mean(+d,2)';
data.mean = pzero;
data.sig  = sig;
data.prec = prec;

%A  = Q * diag(sqrt(abs(L)));      % Data in a pseudo-Euclidean space

% Determine the mapping depending on pzero
if pzero == 0,
  W = mapping(mfilename,'trained',data,[],m,sum(sig));
else
  data.mdis = +d(:,pzero)';
  W = mapping(mfilename,'trained',data,[],m,sum(sig));
end
W = setname(W,'PE embedding');
return
