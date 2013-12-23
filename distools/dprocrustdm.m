% DPROCRUSTDM Distance Matrix between Datasets based on Extended Procrustes Problem
%
%   [D,T] = DPROCRUSTDM(X,Y,SC,EXT)
%
% INPUT
%   X   N x Mx or N x Mx x Kx data
%   Y   N x My or N x My x Ky data; My <= Mx
%   SC  Parameter (1/0) indicating whether to scale the distance to [0,1]
%       by normalizing X and Y by their Frobenius norms (optional; default: 1)
%   EXT Parameter (1/0) indicating extended (1) or orthogonal (0)
%       Procrustes problem (optional; default: 1)
%
% OUTPUT
%   D   Kx x Ky Distance matrix
%   T   Kx x Ky Transformation structure
%
% DEFAULT
%   SC  = 1
%   EXT = 1
%
% DESCRIPTION
% Given two 2D matrices X and Y, extended Procrustes analysis, EXT = 1, finds
% a linear transformation based on shift, orthogonal transformation and scaling
% of the points in Y to fit the points in X.  This is done by minimizing the sum
% of squared differences, which is also the Frobenius norm between X and the
% transformed Yt. Yt = alpha*Y*Q+1*beta^T, where alpha is the scaling scalar,
% Q is the orthogonal transformation, beta is the shift vector and 1 is the
% vector of all ones. If SC = 0, then the resulting difference is returned as
% the dissimilarity D. So, the parameters are found in the least square sense
% such that  min ||X - alpha*Y*Q - 1*beta^T||^2 holds.
% Then, D = norm(X-Yt,'Frobenius').
% If SC = 1, then the resulting distance D is scaled to [0,1] by normalizing it
% by NORM(Xc,'Frobenius'), where Xc is X shifted to the origin.
%
% Orthogonal Procrustes analysis, EXT = 0, neglects alpha and beta and focuses
% on the orthogonal transformation only. So, the above holds for Yt = Y*Q.
%
% X and Y should have the same number of points, as Procrustes analysis matches
% X(i,:) to Y(i,:). If dim(Y) < dim(X), then columns of zeros are added to Y.
% If X and Y are 3D matrices, then we consider sets of 2D matrices to be
% compared, which results in a Kx x Ky distance matrix D.
%
% IMPORTANT
% Note that D = DPROCRUST(X,X,0,1) is assymetric and D = DPROCRUST(X,X,S,EXT)
% is symmetric, otherwise.
%
% T is a structure of the size Kx x Ky with the fields of alpha, beta and Q
% if EXT = 1 and with the field of Q, if EXT = 0. For instance, T(i,j).Q is the
% orthogonal trnasformation of Y(:,:,j) to fit Y(:,:,i).
%
% This routine can be used to match two results of MDS, or KPCA/PSEM for
% an Euclidean embedding, or two shapes (described by contour points) of
% known point correspondences.
%
% SEE ALSO
% MDS, PCA, PSEM, KPCA, MAPPINGS, DATASETS
%
% REFERENCE
% 1. J.C. Gower, Generalized Procrustes analysis. Psychometrika, 40, 33-51, 1975.
% 2. I. Borg, and P. Groenen, Modern Multidimensional Scaling. Springer, New York, 1997.
% 3. http://e-collection.ethbib.ethz.ch/ecol-pool/bericht/bericht_363.pdf

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function [D,T] = dprocrustdm(X,Y,normalize,extPA)

% Let X and Y be first centered at the origin.
% [U,S,V]=eig(X'*Y);
%
% For non-normalized data, the transformation is given as
% Q     = V*U';
% alpha = trace(X'*Y*Q) / trace(Y'*Y) = sum(diag(S)) / norm(Y,'Frobenius')^2
% beta  = (mean(X) - alpha*mean(Y)*Q);

if nargin < 3,
  normalize = 1;
end
if nargin < 4,
  extPA = 1;
end

if normalize~=0 & normalize ~=1,
  error('SC should be either 0 or 1.');
end

if extPA~=0 & extPA~=1,
  error('EXT should be either 0 or 1.');
end

if isdataset(X), X = +X; end
if isdataset(Y), Y = +Y; end

[nx,mx,kx] = size(X);
[ny,my,ky] = size(Y);

if ny ~= nx,
  error('X and Y should have the same number of points.');
end
if my > mx,
  error('Y cannot have more columns (dimensions) than X.');
end
% Add columns of zero if dim(Y) < dim(X)
if my < mx,
  if ky == 1,
    Y = [Y zeros(ny,mx-my)];
  else
    Y(:,my+1:mx,:) = zeros(ny,mx-my,ky);
  end
end

XX = X;
YY = Y;

% Center the data at the origin
Xme = mean(X,1);
Yme = mean(Y,1);
X   = X - repmat(Xme,nx,1);
Y   = Y - repmat(Yme,ny,1);

if normalize,
  % Compute the square Frobenius norm and scale the data
  % For data A, norm(A,'Frob') = sqrt(trace(A'*A)) = sqrt(sum(sum(A.^2))
  Xfn = sqrt(sum(sum(X.^2,1),2));
  Yfn = sqrt(sum(sum(Y.^2,1),2));
  X   = X./repmat(Xfn,nx,mx);
  Y   = Y./repmat(Yfn,ny,my);
  Xfn = squeeze(Xfn);
  Yfn = squeeze(Yfn);
else
  Yfn2 = squeeze(sum(sum(Y.^2,1),2));
end


% Find the optimal transformation parameters of Yt = alpha*YY*Q+1*beta^T,
% rotation matrix Q, scaling alpha and shift vector beta
if normalize,
  for i=1:kx
    for j=1:ky
      G  = X(:,:,i)'*Y(:,:,j);
      [U,S,V] = svd(G);
      Q  = V*U';
      scY(i,j) = sum(diag(S));
      if nargout > 1,
        T(i,j).Q = Q;
        if extPA,
          T(i,j).alpha = scY(i,j) * (Xfn(i) /Yfn(j));
          T(i,j).beta  = squeeze(Xme(:,:,i)) - T(i,j).alpha*squeeze(Yme(:,:,j)) * T(i,j).Q;
        end
      end
    end
  end
  % Distance between X and Yt = alpha*Y*Q+1*beta^T
  D = real(sqrt(1 - scY.^2));
else
  for i=1:kx
    for j=1:ky
      G = X(:,:,i)'*Y(:,:,j);
      [U,S,V] = svd(G);
      Q = V*U';
      scY(i,j) = sum(diag(S));
      if nargout > 1,
        T(i,j).Q = Q;
        if extPA,
          T(i,j).alpha = scY(i,j) ./ Yfn2(j);
          T(i,j).beta  = squeeze(Xme(:,:,i)) - T(i,j).alpha*squeeze(Yme(:,:,j)) * T(i,j).Q;
          D(i,j) = norm(XX(:,:,i) -T(i,j).alpha*YY(:,:,j)*T(i,j).Q+ones(ny,1)*T(i,j).beta,'fro');
        else
          D(i,j) = norm(XX(:,:,i) -YY(:,:,j)*T(i,j).Q,'fro');
        end
      else
        if extPA,
          alpha  = scY(i,j) ./ Yfn2(j);
          beta   = squeeze(Xme(:,:,i)) - alpha*squeeze(Yme(:,:,j)) * Q;
          D(i,j) = norm(XX(:,:,i) - alpha*YY(:,:,j)*Q + ones(ny,1)*beta,'fro');
        else
          D(i,j) = norm(XX(:,:,i) - YY(:,:,j)*Q,'fro');
        end
      end
    end
  end
end

% Take care that distances are real and nonnegative
D = real(D);
D(find (D < 0)) = 0;
if kx == ky & kx > 1 & isymm(D,1e-12),
  D(1:kx+1:end) = 0;
  D = 0.5*(D+D');
end
