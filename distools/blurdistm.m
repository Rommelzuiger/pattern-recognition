%BLURDM Blurred Euclidean distance between blobs
%
%   D = BLURDISTM(A,B,RA,RB)
%
% INPUT
%   A     MYA x MXA x NA matrix of NA binary images of the size MYA x MXA
%   B     MYB x MXB x NB matrix of NA binary images of the size MYB x MXB
%   RA,RB Blurring paramters (optional, default: 1, no blurring)
%
% OUTPUT
%   D     NA x NB Blurred Euclidean distance matrix
%
% DESCRIPTION
% Computes the blurred Euclidean distance matrix D between sets of blobs.
% Images are first uniformly blurred with the sizes of RA x RA or RB x RB.
% Next they are rescaled to square images with the sizes equal to
% max([MYA,MXA,MYB,MXB]) vy using bilinear interpolation. Finally, the
% Euclidean distance is computed.
%
% Preferably use NA <= NB, as it is faster.

% Copyright: R.P.W. Duin, r.duin@ieee.org
% Ela Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function d = blurdistm(A,B,ra,rb,fid)
if nargin < 3,
  ra = 1;
end
if nargin < 4,
  rb = ra;
end
if nargin < 5,
  fid = 0;
end

[ma1,ma2,na] = size(A);
[mb1,mb2,nb] = size(B);
m = max([ma1,ma2,mb1,mb2]);
d = zeros(na,nb);
filta = exp(-[-ra:ra].^2/(0.5*ra*ra));
filta = filta./sum(filta);
filtb = exp(-[-rb:rb].^2/(0.5*rb*rb));
filtb = filtb./sum(filtb);

if fid > 0,
  figure(1); clf;
  figure(2); clf;
end

for i=1:na
  a = A(:,:,i);
  J = find(any(a));
  J = [min(J):max(J)];
  K = find(any(a'));
  K = [min(K):max(K)];
  a = double(a(K,J));

% figure(1); imagesc(a); drawnow

  a = bord(a,0,ra);
  a = conv2(filta,filta,a,'same');
  a = a(ra:end-ra,ra:end-ra);
  a = imresize(a,[32 32],'bilnear');
  if fid > 0,
    figure(1); imagesc(a); drawnow
  end
  for j = 1:nb
    b = B(:,:,j);
    J = find(any(b));
    J = [min(J):max(J)];
    K = find(any(b'));
    K = [min(K):max(K)];
    b = double(b(K,J));
    b = bord(b,0,rb);
    b = conv2(filtb,filtb,b,'same');
    b = b(rb:end-rb,rb:end-rb);
    b = imresize(b,[32 32],'bilnear');
    d(i,j) = sqrt(sum((a(:)-b(:)).^2));
    if fid > 0
      fprintf(fid,'%5d %5d %10.3f \n',i,j,d(i,j));
      figure(2); imagesc(b); drawnow;
    end
  end
end
