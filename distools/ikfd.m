%IKFD Indefinite Kernel Fisher Discriminant
%
%   W = IKFD(K,ALF,BIAS)
%
% INPUT
%   K     NxN kernel or similarity matrix (dataset)
%   ALF   Regularization constant
%   BIAS  Use bias (1) or not (0)
%
% OUTPUT
%   W     Trained Kernel Fisher discriminant
%
% DEFAULT
%   ALF  = 0.0001
%   BIAS = 1
%
% DESCRIPTION
% Finds a Fisher linear discriminant in a kernel-induced space.
% Regularization is necessary.
% Multi-class classifier is trained one-vs-all classes.
%
% SEE ALSO
% KPCA, KFD, FISHERC, MAPPINGS, DATASETS
%
% REFERENCE
% S. Mika, G. Ratsch, J. Weston, B. Scholkopf, and K.-R. Muller.
% Fisher discriminant analysis with kernels. In: Neural Networks for
% Signal Processing IX, pages 41-48, 1999.
%

% Copyright: Ela Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester
%


function W = ikfd(K,alf,isb0)
if nargin < 3,
  isb0 = 1;
end
if nargin < 2 | isempty(alf),
  alf = 0.0001;
end
if nargin < 1 | isempty(K),
  W = mapping(mfilename,{alf,isb0});
  W = setname(W,'IKFD');
  return
end

if alf <= 0,
  error('A small positive regularization ALF is necessary.');
end

islabtype(K,'crisp');
isvaldset(K,1,2);     % At least one object per class and two classes

lab      = getnlab(K);
lablist  = getlablist(K);
[n,k,C]  = getsize(K);


% If more than two classes, train in the one-against-all strategy.
if C > 2,
  W = [];
  for i=1:C
    mlab = 2 - (lab == i);
    KK   = setlabels(K,mlab);
    KK   = remclass(KK,0);
    if ~isempty(K.prior)
      KK = setprior(KK,[K.prior(i),1-K.prior(i)]');
    end

    v  = ikfd(KK,alf,isb0);
    W  = [W,setlabels(v(:,1),lablist(i,:))];
  end
  return

else    % Two classes

%  V  = kcenterm(K);   % Center the kernel
%  KK = +(K*V);        % but ... centering is not necessary

  KK = +K;
  Y  = 3 - 2 * lab;   % Set labels to +/-1
  I1 = find(Y==1);
  I2 = find(Y==-1);
  n1 = length(I1);
  n2 = length(I2);

  M1 = sum(KK(:,I1),2)/n1;
  M2 = sum(KK(:,I2),2)/n2;
  M  = (M1-M2)*(M1-M2)';
  N  = KK(:,I1)*(eye(n1) - ones(n1)/n1)*KK(:,I1)' + KK(:,I2)*(eye(n2) - ones(n2)/n2)*KK(:,I2)';

  % Regularization
  N = N + alf * eye(n);

  % Optimization
  %[W,V,U] = svds(inv(N)*M,1);
  W = inv(N)*(M1-M2);   % The same as above up to scaling

  % Project data on the found direction
  b0 = 0;
  if isb0,
    % Find the free parameter
    b0 = -W'*(M1+M2)/2;
  end
end


% Determine the mapping
W = affine(W,b0,K,lablist,k,2);
%W = cnormc(W,K);
%W = V*W;
W = setname(W,'IKFD');
return;
