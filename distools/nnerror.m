%NNERROR Exact expected NN error from a dissimilarity matrix (2)
%
%   E = NNERROR(D,M)
%
% INPUT
%   D 	NxN dissimilarity dataset
%   M 	Number of objects to be selected per class (optional)
%
% OUTPUT
%   E 	Expected NN errror 
%
% DEFAULT
%   M = minimum class size minus 1
%
% DESCRIPTION
% An exact computation is made of the expected NN error for a random 
% selection of M for training. D should be a dataset containing a labeled 
% square dissimilarity matrix.

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands


function e = nnerror(D,n)
isdataset(D);
discheck(D);
nlab     = getnlab(D);
[m,mm,c] = getsize(D);

% Number of objects per class
nc = classsizes(D);

% Compute for all training set sizes
if nargin < 2
	n = min(nc)-1;
end

if length(n) > 1
	e = zeros(1,length(n));
	for j=1:length(n)
		e(j) = nnerror(D,n(j));
	end
else
	if n >= min(nc)
		error('Requested size of the training set is too large.')
	end
	% Call for the given sample size	
	[D,I] = sort(D);
	I     = reshape(nlab(I),m,m);   % Order objects according to their distances
	ee    = zeros(1,m);
	
	% Loop over all classes
	for j = 1:c
		% Find probabilities Q that objects of other classes are not selected
		Q = ones(m,m);
		for i = 1:c
			if i~=j
				[p,q] = nnprob(nc(i),n); 
				q = [1,q];
				C = cumsum(I==i)+1;
				Q = Q.*reshape(q(C),m,m);
			end
		end
		
		% Find probabilitues P for objects of this class to be the first
		[p,q] = nnprob(nc(j),n); 
		p = [0,p];
		C = cumsum(I==j)+1;
		P = reshape(p(C),m,m);

		% Now estimate the prob EC it is really the NN
		J     = find(I==j);
		EC    = zeros(m,m);
		EC(J) = P(J).*Q(J);
		
		% Determine its error contribution
		L     = find(nlab==j);
		ee(L) = 1-sum(EC(2:m,L))./(1-EC(1,L));	% Correct for the training size
	end

	% Average for the final result
	e = abs(mean(ee));
end



%NNPROB Probability of selection as the nearest neighbor
%
% [P,Q] = NNPROB(M,K)
%
% If K objects are selected out of M, then P(i) is the probability
% that the i-th object is the nearest neigbor and Q(i) is the probability
% that this object is not selected.

function [p,q] = nnprob(m,k)
p = zeros(1,m);
q = zeros(1,m);
q(1) = (m-k)/m;
p(1) = k/m;
for i=2:(m-k+1)
	q(i) = q(i-1)*(m-k-i+1)/(m-i+1);
	p(i) = q(i-1)*k/(m-i+1);
end

