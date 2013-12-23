%KSVC_NU Kernel Support Vector Classifier: NU algorithm
%
%   [W,J,C,REG] = KSVC_NU(K,NU,KC,R)
%
% INPUT
%   K   NxN Kernel dataset
%   NU  Regularization parameter (0 < NU < 1);
%       expected fraction of support vectors (optional; default: 0.25)
%   KC  Kernel centering, 1/0, (optional; default: 1)
%   R   Parameter: -1,0,1,2.
%       -1 or 'flip',   Changes the kernel by flipping negative eigenvalues to
%          positive
%        0 or '',       Uses the kernel as it is
%        1 or 'reg',    Checks positive definiteness and regularizes the kernel by
%          adding the minimal constant in the form of 10^i to the diagonal
%        2 or 'lamreg', Checks positive definiteness and regularizes the kernel by
%          adding the minimal constant to the diagonal (equal to the magnitude of
%          the smallest negative eigenvalue)
%          (optional; default: 0, do not change the kernel)
%
% OUTPUT
%   W   Mapping: Support Vector Classifier
%   J   Object identifiers of support objects
%   C   Equivalent C regularization parameter of the KSVC algorithm
%   REG Regularization parameter added to the diagonal, if used (R=1,2); a vector
%       of eigenvalues of K (R=-1), or -1 if not checked (R=0)
%
% DEFAULT
%   NU  = 1
%   KC  = 1
%   R   = 0 (use the kernel as it is)
%
% DESCRIPTION
% Determines a support vector machine (SVM) for the kernel K. Quadratic programming
% formulation is used to solve the problem. K can be positive definite or indefinite.
% indefinite. J is a list of the indices of the support objects from K. NU belongs
% to the interval (0,1). The larger NU, the larger class overlap. Default NU = 0.25.
%
% For positive semidefinite kernels, NU is bounded from above by NU_MAX, where
% NU_MAX = (1 - ABS(Lp-Lm)/(Lp+Lm)), where Lp (Lm) is the number of positive
% (negative) samples. If NU > NU_MAX is supplied to the routine it will be changed to
% the NU_MAX.
% If NU is less than some NU_MIN which depends on the overlap between the classes,
% the algorithm will typically take long time to converge (if at all). So, it is
% advisable to set NU larger than the expected overlap.
% The output is rescaled in such a way as if returned by KSVC with the parameter C.
%
% If R = 2, then K is regularized by adding the smallest constant possible to the
% diagonal to make it positive definite. If R = 1, then K is regularized by adding
% the smallest constant in the form of 10^i to the diagonal. If R = -1, then the
% eigendecomposition of K is found as K = Q*L*Q' and the negative eigenvalues are
% flipped, so the SVM is built on K = Q*abs(L)*Q'.
%
% IMPORTANT
% The classifier cannot always be constructed for an indefinite kernel. If the norm
% of the determined weight vector V in the Krein (pseudo-Euclidean) space induced
% by K is negative, i.e. if V'*K*V < 0, then the proper SVM cannot be built. This
% happens when the hyperplane in the Krein space lies 'more' in the negative subspace
% than in the positive one. This means that the kernel is strongly indefinite.
% Currently, a warning is given and a pseudo-Fisher classifier is trained, instead.
% The pseudo-Fisher classifier is also trained in any situation when the quadratic
% optimization fails to find a solution.
%
% REMARKS
% Note that if D is a symmetric distance/dissimilarity matrix, then K = -D is an
% (indefinite) kernel. If D.^2 is a square Euclidean distance matrix, then K = -D.^2
% is a proper (conditionally negative definite) kernel. So, a linear SVM on some
% data X, based on the kernel K = X*X' is the same classifier as the kernel-based
% SVM on K = -D(X,X).^2.

% DESCRIPTION
% Determines a support vector machine (SVM) for the kernel K. Quadratic programming
% formulation is used
% The smaller C, e.g.
% C < 1, the larger class overlap imposed. Optimal C will be different for a regularized
% kernel and for an indefinite one.
%
% SEE ALSO
% KSVO_NU, KSVO, KSVC, MAPPINGS
%
% LITERATURE
% 1. B.Scholkopf, A. Smola, Learning with kernels, MIT press, 2001,
%    http://www.learning-with-kernels.org/.
% 2. B. Haasdonk, Feature Space Interpretation of SVMs with Indefinite Kernels.
%    IEEE Trans. on PAMI, 27(4):482-492, 2005.
% 3. E. Pekalska, P. Paclik, R.P.W. Duin,  A Generalized Kernel Approach to
%    Dissimilarity-based Classification, JMLR, vol.2, no.2, 175-211, 2002.


% Copyright: Elzbieta Pekalska, Robert P.W. Duin, ela.pekalska@googlemail.com
% Based on SVC.M by D.M.J. Tax, D. de Ridder and R.P.W. Duin and SVC_NU.M by S.Verzakov
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [W,J,reg,err,K] = ksvc(K,nu,kc,r,no)
prtrace(mfilename);

if nargin < 5,
  % Multi-class problems are solved by one-vs-all by calling KSVC.
  % no is the class number in such a situation; 0 is the standard 2-class case
  no = 0;
end
if nargin < 4,
  r = 0;
end
if nargin < 3,
  kc = 1;
end
if nargin < 2 | isempty(nu),
  nu = 0.25;
  prwarning(3,'Regularization parameter NU set to 0.25\n');
end
if nargin < 1 | isempty(K),
  W = mapping(mfilename,{nu,kc,r});
  W = setname(W,'Kernel Support Vector Classifier (KSVC-NU)');
  return;
end


if all(kc ~= [0 1]),
  error('Wrong KC parameter.');
end

switch lower(r)
  case {'flip', -1},  r = -1;
  case {'', 0},       r = 0;
  case {'reg', 1},    r = 1;
  case {'lamreg', 2}, r = 2;
  otherwise
    error('Wrong parameter R.');
end


% TRAIN THE CLASSIFIER
if ~isa(nu,'mapping')
  islabtype(K,'crisp');
  isvaldset(K,1,2);      % Expect at least 1 object per class and 2 classes
  [m,k,c] = getsize(K);
  if m ~=k,
    error('The kernel is not a square matrix.');
  end
  nlab = getnlab(K);

  % The SVC is basically a two-class classifier.
  % Multi-class problems are trained one-versus-all.

  if c == 2   % Two-class classifier
  % Compute the parameters for the optimization:
    y = 3 - 2*nlab;    % y = +1/-1

    prec = 1e-12;
    if ~issym(K,prec),
      prwarning(1, 'The kernel is not symmetric. The values are averaged out.')
      K = (K+K')/2;
    end

    if kc,   % Center the data
      me = mean(+K,2)';                % Store the mean value

      B = -repmat(1/m,m,m);
      B(1:m+1:end) = B(1:m+1:end)+1;   % B = eye(m) - ones(m,m)/m
      K = B*K*B;                       % K is now centered
    else
      me = [];
    end
    K = (K+K')/2;                      % K is averaged out as QP solver is sensitive to small asymmetry


    % Check feasibility of the kernel:
    if r == 0 & dmeans(K,y) < 0,
      if no > 0,
        ss = [num2str(no) '-vs-all '];
      else,
        ss = '';
      end
      prwarning(1,['The kernel is badly indefinite. The ' ss ' SVM cannot be defined. Pseudo-Fisher is computed instead.']);
      err = 1;
      v   = pinv([K ones(m,1)])*y;
      J   = [1:m]';
      T   = [];
      reg = [];
    else
      % Perform the optimization
      [v,J,T,reg,err] = ksvo_nu(+K,y,nu,r);
    end

    % Store the results
    W = mapping(mfilename,'trained',{me,J,T,v},getlablist(K),k,2);
    %W = cnormc(W,a);
    W = setname(W,'Kernel Support Vector Classifier (KSVC-NU)');
    W = setcost(W,K);
    % J = K.ident(J);

  else   % MULTI-CLASS CLASSIFIER: C > 2

    % MCLASSC cannot be used here as we have a kernel K
    W = [];
    J = zeros(m,1);
    lablist = getlablist(K);

    for i=1:c
      lab = 2 - (nlab == i);   % lab = 1/2
      KK  = setlabels(K,lab);
      KK  = remclass(KK,0);
      KK  = setfeatlab(K,lab);
      if ~isempty(K.prior)
        KK = setprior(KK,[K.prior(i),1-K.prior(i)]');
      end
      [V,j,reg(i),err(i)]= ksvc_nu(KK,nu,kc,r,i);
      W = [W,setlabels(V(:,1),lablist(i,:))];
      J(j) = 1;
    end
    J = find(J);
  end

else

% EXECUTE THE CLASSIFIER
  % nu is an SVM classifier now
  n = size(K,1);
  w = +nu;

  % Get the parameters from the classifier
  me = w{1};
  J  = w{2};

  % The first   % The first parameter w{1} stores the mean of the kernel.
  % When it is non-empty, data centering should also be applied
  % to the test kernel.

  if ~isempty(me),
    % Center the kernel
    m = length(me);
    B = -repmat(1/m,m,m);
    B(1:m+1:end) = B(1:m+1:end) + 1;    % B = eye(m) - ones(m,m)/m
    K = (K - me(ones(n,1),:)) * B;
  end

  if ~isempty(w{3}),   % this is the transformation that reverses the negative eigenvalues
    K = K*w{3};
  end

  % The classifier is stored in w{4}
  % Data is mapped by the kernel, now we just have a linear classifier  w*x+b:
  d = [+K(:,J) ones(n,1)] * w{4};
  d = sigm([d -d]);
  W = setdat(K,d,nu);
end

return




function dm = dmeans(K,y)
% Computes the square pseudo-Euclidean distance between the
% means of the two classes. This can be done by using the kernel only.
% Negative value means that the contribution of the negative
% subspace in the pseudo-Euclidean sense is too large.

yy = y;
Z  = find (y == 1);
T  = find (y == -1);
yy(Z) = 1/length(Z);
yy(T) = -1/length(T);
dm    = yy'*(+K)*yy;
return;
