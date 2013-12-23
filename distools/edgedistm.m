%EDGEDISTM  Distance Matrix between Images based on Their Edges
%
%   D = EDGEDISTM (A,B,DIST,K)
%
% INPUT
%   A     Matrix of N binary images
%   B     Matrix of M binary images
%   DIST  Distance measure (optional; default: 'MHAUS')
%         'HAUS'  - Hausdorff distance
%         'MHAUS' - modified Hausdorff distance
%         'R75'   - based on the 75th quantile
%         'R90'   - based on the 90th quantile
%   K     Number of rows in an image if images are stored as row vectors
%         (optional; default: sqrt(size(A,1)))
%
% OUTPUT
%   D     NxM dissimilarity matrix or dataset
%
% DESCRIPTION
% Computes the distance matrix D between two sets of binary images, A and B.
% First edges are extracted, then a distance defined by DIST is computed.
% If the size of A is N x P and the size of B is M x P, then images are
% stored as row vectors. P should be such that P = K*L, where K and L are
% image sizes.
%
% If the size of A is Kx x LX x N and the size of B is KY x LY x M, then
% images are stored as 2D arrays.
%
% If A and B are datasets, then D is a dataset as well with the labels defined
% by the labels of A and the feature labels defined by the labels of B. If A is
% not a dataset, but a matrix of doubles, then D is also a matrix of doubles.

% Copyright: Robert P.W. Duin, r.p.w.duin@prtoos.org, Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function D = edgedistm (A,B,dist,k)

if nargin < 3,
  dist = 'MHAUS';
else
  dist = upper(dist);
end

if nargin < 2,
  B = A;
end

switch ndims(A)
  case 2,
    isda = isdataset(A);
    isdb = isdataset(B);
    a = +A;
    b = +B;
    [ma, na] = size(a);

    if nargin < 4,
      k = round(sqrt(na));
    end
    l = round(na/k);

    if (k * l) ~= na,
      error ('Sizes of the images are incorrect.');
    end

    [mb, nb] = size(b);
    if na ~= nb,
      error ('Number of features should be the same.');
    end

    D = zeros(ma,mb);
    for i=1:ma
      for j=1:mb
        [Ia,Ja] = find (edge (reshape (a(i,:),k,l),'canny'));
        [Ib,Jb] = find (edge (reshape (b(j,:),k,l),'canny'));
        D(i,j)  = edist([Ia Ja],[Ib Jb],dist);
      end
    end

  % Set object labels and feature labels
  if xor(isda, isdb),
    prwarning(1,'One matrix is a dataset and the other not. ')
  end
  if isda,
    if isdb,
      D = setdata(A,D,getlab(B));
    else
      D = setdata(A,D);
    end
    D.name = 'Distance matrix';
    if ~isempty(A.name)
      D.name = [D.name ' for ' A.name];
    end
  end

  case 3,
    [ka,la,ma] = size(a);
    [kb,lb,mb] = size(b);

    D = zeros(ma,mb);
    for i=1:ma
      for j=1:mb
        [Ia,Ja] = find (edge (a(:,:,i),'canny'));
        [Ib,Jb] = find (edge (b(:,:,j),'canny'));
        D(i,j)  = edist ([Ia Ja],[Ib Jb],dist);
      end
    end

  otherwise
    error('The number of matrix''s dimensions should be 2 or 3.');
end

D(find(D < eps)) = 0;
return








function dd = edist (aa,bb,dist)
dab   = sqrt(distm(aa,bb));
minab = min(dab);

dba   = sqrt(distm(bb,aa));
minba = min(dba);

switch dist,
  case 'haus',
    dAB = max(minab);
    dBA = max(minba);

  case 'mhaus'
    dAB = mean(minab);
    dBA = mean(minba);

  case {'r75','r90'}
    if strcmp(dist,'r75')
      r = 75;
    else
      r = 90;
    end;
    minab = sort(minab);
    dAB   = minab(round (r/ma));
    minba = sort(minba);
    dBA   = minba(round (r/mb));
  otherwise
    error ('Wrong distance measure.')
end

dd = max(dBA, dAB);
return
