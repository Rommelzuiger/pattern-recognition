%KSVC Kernel Support Vector Classifier on a kernel matrix
%
%   [W,J,REG] = KSVC(K,C,KC,R)
%
% INPUT
%   K   NxN Kernel dataset
%   C   Regularization parameter (optional; default: 1)
%   KC  Kernel centering, 1/0, (optional; default: 1)
%   R   Parameter: -1,0,1,2
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
%   REG Regularization parameter added to the diagonal, if used (R=1,2); a vector
%       of eigenvalues of K (R=-1), or -1 if not checked (R=0)
%
% DEFAULT
%   C  = 1
%   KC = 1
%   R  = 0 (use the kernel as it is)
%
% DESCRIPTION
% Determines a support vector machine (SVM) for the kernel K. Quadratic programming
% formulation is used to solve the problem. K can be positive definite or indefinite.
% J is a list of the indices of the support objects from K. The smaller C, e.g.
% C < 1, the larger class overlap imposed. Optimal C will be different for regularized
% and indefinite kernels.
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
% happens when the hyperplane in the Krein space lies 'more' in the negative
% subspace than in the positive one. This means that the kernel is strongly
% indefinite. Currently, a warning is given and a pseudo-Fisher classifier is trained,
% instead. The pseudo-Fisher classifier is also trained in any situation when the
% quadratic optimization fails to find a solution.
%
% REMARK
% Note that if D is a symmetric distance/dissimilarity matrix, then K = -D is an
% (indefinite) kernel. If D.^2 is a square Euclidean distance matrix, then K = -D.^2
% is a proper (conditionally negative definite) kernel. So, a linear SVM on the
% vectorial data X, based on the kernel K = X*X', is equivalent to the kernel-based
% SVM on K = -D(X,X).^2.
%
% SEE ALSO
% MAPPINGS, DATASETS, KSVO, SVC,
%
% LITERATURE
% 1. B.Scholkopf, A. Smola, Learning with kernels, MIT press, 2001,
%    http://www.learning-with-kernels.org/.
% 2. B. Haasdonk, Feature Space Interpretation of SVMs with Indefinite Kernels.
%    IEEE Trans. on PAMI, 27(4):482-492, 2005.
% 3. E. Pekalska, P. Paclik, R.P.W. Duin,  A Generalized Kernel Approach to
%    Dissimilarity-based Classification, JMLR, vol.2, no.2, 175-211, 2002.

% Copyright: Elzbieta Pekalska, Robert P.W. Duin, ela.pekalska@googlemail.com
% Based on SVC.M by D.M.J. Tax, D. de Ridder and R.P.W. Duin
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [W,J,reg,err,K] = ksvc(K,C,kc,r,no)
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
if nargin < 2 | isempty(C),
  C = 1;
  prwarning(3,'Regularization parameter C set to 1.\n');
end
if nargin < 1 | isempty(K),
  W = mapping(mfilename,{C,kc,r});
  W = setname(W,'Kernel Support Vector Classifier (KSVC)');
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
if ~isa(C,'mapping')
  islabtype(K,'crisp');
  isvaldset(K,1,2);      % Expect at least 1 object per class and 2 classes
  [m,k,c] = getsize(K);
  if m ~=k,
    error('The kernel is not a square matrix.');
  end
  nlab  = getnlab(K);

  % The SVM is basically a two-class classifier.
  % Multi-class problems are trained one-versus-all.

  if c == 2   % Two-class classifier
    % Compute the parameters for the optimization:
    y = 3 - 2*nlab;

    prec = 1e-12;
    if ~issym(K,prec),
      prwarning(1, 'The kernel is not symmetric. The values are averaged out.')
    end

    if kc,   % Center the data
      me = mean(+K,2)';                % Store the mean value

      B = -repmat(1/m,m,m);
      B(1:m+1:end) = B(1:m+1:end)+1;   % B = eye(m) - ones(m,m)/m
      K = B*K*B;                       % K is now centered
    else
      me = [];
    end
    K = (K+K')/2;                      % K is averaged as QP solver is sensitive to small asymmetry


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
      [v,J,T,reg,err] = ksvo(+K,y,C,r);
    end

    % Store the results
    W = mapping(mfilename,'trained',{me,J,T,v,reg},getlablist(K),k,2);
    % W = cnormc(W,K);
    W = setname(W,'Kernel Support Vector Classifier (KSVC)');
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
      [V,j,reg(i),err(i)]= ksvc(KK,C,kc,r,i);
      W = [W,setlabels(V(:,1),lablist(i,:))];
      J(j) = 1;
    end
    J = find(J);
  end

else

% EXECUTE THE CLASSIFIER
  % C is a SVM classifier now
  n = size(K,1);
  w = +C;

  % Get the parameters from the classifier
  me = w{1};
  J  = w{2};

  % The first parameter w{1} stores the mean of the kernel.
  % When it is non-empty, kernel centering should also be applied
  % to the test kernel.

  if ~isempty(me),
    m = length(me);
    % Center the kernel
    B = -repmat(1/m,m,m);
    B(1:m+1:end) = B(1:m+1:end) + 1;    % B = eye(m) - ones(m,m)/m
    K = (K - me(ones(n,1),:)) * B;
  end

  if ~isempty(w{3}),   % this is the transformation that reverses the negative eigenvalues
    K = K*w{3};
  end

  % The classifier is stored in w{4}
  % Data is mapped by the kernel, now we just have a linear classifier  w*x+b:
  d = [K(:,J) ones(n,1)] * w{4};
  d = sigm([d -d]);
  W = setdat(K,d,C);
end
return




function dm = dmeans(K,y)
% Coputes the square pseudo Euclidean distance between the
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
