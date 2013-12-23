%BLOBBOX Find box around a binary blob and resample
%
%   B = BLOBBOX(A,ny,nx)
%
% INPUT
%   A     MY x MX x N matrix of N binary images of the size MY x MX
%   NX,NY Resampling sizes (optional, default: 16)
%
% OUTPUT
%   B     NY x NX x N matrix of N binary images of the size NY x NX
%
% DEFAULT
% NX = NY = 16.
%
% DESCRIPTION
% For an MY x MX x N set of N binary images A the bounding boxes around
% the blobs are computed and resampled with NX x NY pixels.
%

% Copyright: R.P.W. Duin, r.duin@ieee.org
% and Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function b = blobbox(a,ny,nx);
if nargin < 2,
  ny = 16;
end
if nargin < 3,
  nx = ny;
end

[mx,my,n] = size(a);
b = zeros(nx,ny,n);

for i=1:n
  c = a(:,:,i);
  J = find(any(c));
  J = [min(J):max(J)];
  K = find(any(c'));
  K = [min(K):max(K)];
  c = double(c(K,J));
  if length(c(:)) > 0,
    c = bord(c,0);
    b(:,:,i) = imresize(c,[ny,nx]);
  end
end



function C = bord(A,n,m);
% C = bord(A,n,m)
% Puts a border of the width m (default m=1) around the image A
% and gives it value n. If n = NaN: mirror image values.

if nargin == 2;
  m=1;
end
[x,y] = size(A);
if m > min(x,y)
  mm = min(x,y);
  C = bord(A,n,mm);
  C = bord(C,n,m-mm);
  return
end

if isnan(n)
  C = [A(:,m:-1:1),A,A(:,y:-1:y-m+1)];
  C = [C(m:-1:1,:);C;C(x:-1:x-m+1,:)];
else
  bx = ones(x,m)*n;
  by = ones(m,y+2*m)*n;
  C = [by;[bx,A,bx];by];
end
