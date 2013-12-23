%CUTEIGS Select significant eigenvalues from a list
%
%   J = CUTEIGS(L,PREC)
%
% INPUT
%   L     List of eigenvalues
%   PREC  Precision parameter (optional; default: 0.0001)
%
% OUTPUT
%   J   Index of selected eigenvalues
%
% DESCRIPTION
% This is a low-level routine for SELEIGS, which serves PSEM, KPSEM, KPCA and
% PEPCA. It makes use of PCHIP determing piecewise cubic Hermite interpolating
% polynomial P, built on [1:length(LL) LL], where LL = -sort(-abs(L).^3) folloowed
% by a normalization to [0,1]. A third power is used to emphasize the differences
% better, as sometimes there is a very long tail of slowly dropping eigenvalues.
% The number M of significant eigenvalues is determined such that the ratio
% of the estimated derivative of P versus the largest derivative value is smaller
% then the given precision, PREC * cumsum(LL).
%
% Note that M will differ between different samples and for different sizes of L.
%
% SEE ALSO
% MAPPINGS, DATASETS, PSEM, KPSEM, KPCA, PEPCA
%

% Elzbieta Pekalska, e.pekalska@ewi.tudelft.nl
% Faculty of Electrical Engineering, Mathematics and Computer Science,
% Delft University of Technology, The Netherlands


function [J,prec] = cuteigs(L,prec)

if nargin < 2 | isempty(prec),
  prec = 0.0001;
end

tol = 1e-14;

n  = length(L);
LL = L.^3;
[lambda,J] = sort(-abs(LL)/sum(abs(LL)));
lambda = -lambda;

if n > 1,
  M = findJ(lambda,prec);
else
  M = n;
end
JJ = J(1:M);

I1 = find (L(JJ) >= 0);
I2 = find (L(JJ) < 0);

J  = [JJ(I1); JJ(I2)];







function M = findJ (lambda,prec)

n = length(lambda);
if n >= 400,
  Z = [1:80 81:4:n];
elseif n >= 200,
  Z = [1:50 51:2:n];
else
  Z =1:n;
end

if n > 3
  % Hermite piecewise cubic interpolating polynomial
  H = pchip(Z,lambda(Z));

  % Derivative of the Hermite piecewise cubic polynomial
  Hder = H;
  Hder.order = 3;
  Hder.coefs = zeros(H.pieces,3);
  for i=1:Hder.pieces
    pp = polyder(H.coefs(i,:));
    Hder.coefs(i,4-length(pp):3) = pp;
  end
  pder= -ppval(1:n,Hder);   % Evaluated derivative
else
  pder= [lambda(1) -diff(lambda')];   % Approximated derivative
end

% pder(1) is an extrapolated value; pder(2) is used
I  = find(pder/pder(2) < prec*cumsum(lambda'));
II = diff(I);
Q  = find(II > 1);
if ~isempty(Q),
  I = I(Q(1));
end
if ~isempty(I),
  M = I(1)-1;
else
  if n <= 25,
    M = n;
  else
    M = n-1;
  end
end
if M < 1, M = 1; end
