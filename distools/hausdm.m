%HAUSDM Hausdorff and modified Hausdorff distance between datasets of image blobs
%
%  [DH,DM] = HAUSDM(A,B,FID)
%
% INPUT
%   A     XAxYAxNA matrix of NA binary images of the size XA x YA
%   B     XBxYBxNB matrix of NB binary images of the size XB x YB
%   FID   0/1 Report progress on the screen (default: 0)
%
% OUTPUT
%   DH    NAxNB Hausdorff distance matrix
%   DM    NAxNB Modified Hausdorff distance matrix
%
% DESCRIPTION
% Computes a Hausdorff distance matrix DH and a modified Hausdorff distance
% matrix DM between the sets of binary images A and B, or datasets containing 
% them as features. Preferably, NA <= NB (faster computation).
% Progress is reported in fid (fid = 1: on the sreeen).
%
% LITERATURE
% M.-P. Dubuisson and A.K. Jain, "Modified Hausdorff distance for object matching", 
% International Conference on Pattern Recognition, vol. 1, 566-568, 1994.
%

% Copyright: R.P.W. Duin, r.duin@ieee.org
% Faculty of EWI, Delft University of Technology


function [dh,dm] = hausdm(A,B,fid)

if nargin < 3, fid = 0; end
if isdataset(A) & isdataset(B)
	[dh,dm] = hausdm(data2im(A),data2im(B),fid);
	dh = setdata(A,dh);
	dm = setdata(A,dm);
	return
end

[ma1,ma2,na] = size(A);
[mb1,mb2,nb] = size(B);
dh = zeros(na,nb);
dm = zeros(na,nb);
for i=1:na
	a = A(:,:,i); 
	J = find(any(a));
	J = [min(J):max(J)];
	K = find(any(a'));
	K = [min(K):max(K)];
	a = double(a(K,J));
	if length(a(:)) > 0
		a = bord(a,0);
	end
	ca = contourc(a,[0.5,0.5]);
	J = find(ca(1,:) == 0.5);
	ca(:,[J J+1]) =[];
	ca = ca - repmat([1.5;1.5],1,size(ca,2));
	ca = ca/max(ca(:));
	ca = ca - repmat(max(ca,[],2)/2,1,size(ca,2));
	for j = 1:nb
		b = B(:,:,j);
		J = find(any(b));
		J = [min(J):max(J)];
		K = find(any(b'));
		K = [min(K):max(K)];
		b = double(b(K,J));
		if length(b(:)) > 0
			b = bord(b,0);
		end
		cb = contourc(b,[0.5,0.5]);
		J = find(cb(1,:) == 0.5);
		cb(:,[J J+1]) =[];
		cb = cb - repmat([1.5;1.5],1,size(cb,2));
		cb = cb/max(cb(:));
		cb = cb - repmat(max(cb,[],2)/2,1,size(cb,2));
		dab = sqrt(distm(ca',cb'));
		dh(i,j) = max(max(min(dab)),max(min(dab')));
		dm(i,j) = max(mean(min(dab)),mean(min(dab')));
		if fid, disp([i,j,dh(i,j),dm(i,j)]); end
%		fprintf(fid,'%5d %5d %10.3f %8.3f \n',i,j,dh(i,j),dm(i,j));
	end
end
	
