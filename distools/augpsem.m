%AUGPSEM Augmented Pseudo-Euclidean Linear Embedding
%
% --- THIS IS A TEST VERSION! DO NOT RELY ON IT. ----
%
%   [W,SIG,L] = AUGPSEM(D,ALF,P)
%
% INPUT
%   D   NxN symmetric dissimilarity matrix (dataset)
%   ALF Parameter determining the dimensionality and the mapping (optional, defaulf: Inf)
%       (0,1)   - fraction of the total (absolute value) preserved variance
%        Inf    - no dimensionality reduction, keeping all dimensions (it's noisy)
%       'p'     - projection into a Euclidean space based on positive eigenvalues only
%       'PARp'  - projection into a Euclidean space based on the PAR fraction of
%                 positive eigenvalues; e.g. ALF = '0.9p'
%       'n'     - projection into a Euclidean space based on negative eigenvalues only
%       'PARn'  - projection into a Euclidean space based on the PAR fraction of
%                 positive eigenvalues; e.g. ALF = '0.7n'
%       1 .. N  - number of dimensions <= N
%   P   Integer between 0 and N specifying which object is mapped at the origin;
%       0 stands for the mean; (optional, default: 0)
%
% OUTPUT
%   W   Augmented pseudo-Euclidean embedding
%   SIG Signature of the space
%   L   List of eigenvalues
%
% DEFAULT
%   ALF = INF
%   P   = 0
%
% DESCRIPTION
% Linear mapping W onto an M-dimensional pseudo-Euclidean subspace from a symmetric,
% square dissimilarity matrix D such that the dissimilarities are preserved.
% M is determined by ALF. E.g., the subspace is found such that at least a fraction
% ALF of the total variance is preserved for ALF in (0,1). The resulting X is found
% by D*W.
% This is an augmented embedding of the standard pseudo-Euclidean embedding (PSEM)
% such that the dissimilairites are perfectly preserved. The augmentation is by one
% (Euclidean case) or two (pseudo-Euclidean case) dimensions. This has no effect
% on the original data D(R,R); the additional dimension(s) are zeros. The effect
% is apparent when test data D(Te,R) are projected. Due to the projection error,
% the new dissimilarities are not preserved by PSEM, while they are preserved by
% the augmented embedding.
%
% The parameter SIG describes the signature of the subspace. L is a sorted list
% of eigenvalues describing the variances in the (pseudo-)Euclidean space.
%
% SEE ALSO
% PSEM, MAPPINGS, DATASETS
%
% LITERATURE
% 1. A. Harol, E.Pekalska, S.Verzakov, R.P.W. Duin, "Augmented embedding of
% dissimilarity data into (pseudo-)Euclidean spaces", Joint Workshops on S+SSPR
% Lecture Notes in Computer Science, vol. 4109, 613-621, 2006.
% 2. L. Goldfarb, A unified approach to pattern recognition, Pattern Recognition,
% vol.17, 575-582, 1984.

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function [W,sig,L] = augpsem(d,alf,pzero)

if nargin < 3 | isempty(pzero), pzero = 0; end
if nargin < 2 | isempty(alf), alf = inf; end
if nargin < 1 | isempty(d)
   W = mapping(mfilename,alf,pzero);
   W = setname(W,'Augmented PE embedding');
   return
end



if (isdataset(d) | isa(d,'double'))
  if ismapping(alf)

    % APPLY THE MAPPING
    pars = getdata(alf);
    [m,n] = size(d);
    ds    = sum(diag(+d));

    Wp    = pars{1};   % projection onto positive Euclidean subspace
    Wn    = pars{2};   % projection onto negative Euclidean subspace
    dme   = pars{3};   % average square distance of the training disssimilarities
    pzero = pars{4};   % pzero=0 -> the mean of the embedded configuration lies at 0, otherwise the mean lies at pzero
    tsig  = pars{5};   % true signature of PE embedding
    sig   = pars{6};   % signature of the augmented PE embedding
    tol   = pars{7};

    Xp = d*Wp;
    if tsig(2) > 0,
      Xn = d*Wn;
      J  = diag([ones(tsig(1),1); -ones(tsig(2),1)]);
    else
      Xn = [];
      J = diag(ones(tsig(1),1));
    end
    XX = [+Xp +Xn];

    % Compute the true PE norm depending on pzero
    if pzero > 0,
      Xtnorm   = +d(:,pzero).^2;
    else
      Xtnorm   = mean(+d.^2,2) - 0.5*dme;
    end

    % Compute the PE norm of the projected data
    Xpnorm = diag(XX*J*XX');

    % Compute the square projection error
    % In a perfect Euclidean embedding, the true norm should >= the projected norm.
    perr   = (Xtnorm - Xpnorm);  % negative signs of perr indicate PE space

    % To remove noise, set very small values to zero
    P = find(abs(perr) < tol);
    perr(P) = 0;

    % Additional 'positive' dimension is extablish by positive projection error
    Ip =  find(perr >= 0);
    Zp = zeros(m,1);
    Zp(Ip,1) = sqrt(perr(Ip));

    W = Xp;
    if sig(2) > 0,
      % Additional 'negative' dimension is extablish by negative projection error
      In =  find(perr < 0);
      Zn = zeros(m,1);
      Zn(In,1) = -sqrt(abs(perr(In)));
      W = setdat(W, [Xp Zp Xn Zn]);
      if max(abs([Zp; Zn])) < tol & ds > tol & ~issym(d),
         prwarning(1,'Augmented dimensions are close to zero. Perform PSEM, instead.');
      end
    else
      W = setdat(W, [Xp Zp]);
      if max(abs(Zp)) < tol & ds > tol & ~issym(d),
         prwarning(1,'Augmented dimension is close to zero. Perform PSEM, instead.');
      end
    end
    W.user = sig;
    return
  end
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

[W,tsig,L,Q] = psem(d,alf,pzero);
Wp  = psem(W,'p');
Wn  = [];
sig = tsig;   % this is the signature describing augmented mapping
I   = tsig > 0;
sig (I) = sig(I) + 1;
if tsig(2) > 0,
  Wn  = psem(W,'n');
end
dme = mean(+d(:).^2);


%%%%% ------------ THIS IS NOT USED now ------------ %%%%%
%% Compute the true PE norm depending on pzero
%if pzero > 0,
%  Xtnorm = +d(:,pzero).^2;
%else
%  Xtnorm = mean(+d.^2,2) - 0.5*dme;
%end
%
%% Compute the PE norm of the projected data
%% Note that X = Q * diag(sqrt(abs(L))), so
%% Xpnorm = diag (X*J*X') = diag(Q * diag(L) * Q')
%Xpnorm = diag(Q * diag(L) * Q');
%
%% Compute the square projection error
%% It is different from zero if ALF is not 0.99999999 or Inf
%perr   = Xtnorm - Xpnorm;  % negative values indicate PE space
%if sum(perr) < tol,
%  perr = zeros(size(d,1),1);
%end
%%%%% ---------------------------------------------- %%%%%

W = mapping(mfilename,'trained',{Wp,Wn,dme,pzero,tsig,sig,tol},[],m,sum(sig));
W = setname(W,'Augmented PE embedding');
return
