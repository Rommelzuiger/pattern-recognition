%DISRC Dissimilarity based classifier by distance ranking
%
%	W = disrc(A)
%	W = A*disrc
%
% A should be a dissimilarity based dataset with class labels as
% feature labels. The classifier W is based on distance ranking.
% A new dataset B defined by dissimilarities to the same representation
% set A is classified by D = B*W by comparing the dissimilarities in B
% with the ranked dissimilarities in A for each object in the
% representation set separately. The columns in D yield for each
% object in the representation set a distance such that D*classc
% generates for these objects the ranked based probabilities.
%
% Formally W is a set of parallel mappings, so B*W*classd combines
% them by assigning objects to the class with the highest probability.
% Other combiners may be used as well, e.g. B*W*classc*prodc,
% combining all probabilities by the product rule.
%
% See mappings, classc, classd.
%
% $Id: disrc.m,v 1.2 2001/07/10 16:35:58 pavel Exp $

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function w = disrc(a,v)

if nargin < 1 | isempty(a) % empty call (untrained classifier)
        w = mapping('disrc');
		
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
	for j = 1:k
		b = {+a(find(nlab==nf(j)),j) +a(find(nlab~=nf(j)),j)};
		w = [w; mapping('disrc',b,fl(nf(j),:),1,1,1)];
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

