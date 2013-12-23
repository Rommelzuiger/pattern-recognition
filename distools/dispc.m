%DISPC Dissimilarity based Parzen classifier
%
%	W = dispc(A)
%
% A should be a dissimilarity based dataset with class labels as
% feature labels.
%
% $Id: dispc.m,v 1.2 2001/07/10 16:35:58 pavel Exp $

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function w = dispc(a,v)

if nargin < 1 | isempty(a) % empty call (untrained classifier)
        w = mapping('dispc');
		
elseif nargin == 1    % training
	[nlab,lablist,m,k,c,p] = dataset(a);
	[nn,nf,fl] = renumlab(lablist,getfeat(a));
	if c < 2
		error('Dataset should containf more than one class')
	end
	if max(nf) > c
		error('Feature labels of dataset do not match with class labels')
	end
	w = [];
	for j=1:c
		F = find(nf==j);
		L = find(nlab==j);
		n = ceil(sqrt(length(L)));
		[b,r] = sort(+a(L,F));
		for f=1:length(F)
			h(F(f)) = b(r(n,f),F(f));
		end	

		b = {+a(find(nlab==nf(j)),j) +a(find(nlab~=nf(j)),j)};
		w = [w; mapping('distrc',b,fl(nf(j),:),1,1,1)];
	end

%	w = cnormc(w,a); forget normalisation and apriori probs for the moment

elseif nargin == 2 	% evaluation
	[m,k] = size(a);
	b = +v; b1 = b{1}; b2 = b{2};
	n1 = length(b1); n2 = length(b2);
	w1 = zeros(m,1); w2 = zeros(m,1);
	[d R1] = sort([b1;+a]);
	[d R2] = sort([b2;+a]);
	L1 = find(R1>n1);
	L2 = find(R2>n2);
	w1(R1(L1)-n1) = (n1+2-(L1-[1:m]'))/(n1+2);
	w2(R2(L2)-n2) = (n2+2-(L2-[1:m]'))/(n2+2);
%	disp([w1;w2])
	w = invsig(dataset((w1+w2)/2,getlab(a),getfeat(v)));
else
	error('Illegal call')
end

