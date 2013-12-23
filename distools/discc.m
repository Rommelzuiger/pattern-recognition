%DISCC Dissimilarity based combining classifier
%
%	W = discc(D,par)
%
% If D is a square dissimilarity matrix then the distances to each
% representation object are used for building a classifier such that
% sigmoid(a*d+b), can be used to estimate the posterior probability of a
% new object at distance d. The scalar parameters a and b are optimized
% over the training set using a maximum likelihood estimator.
% W is a set of parallel classifiers such that D*W*classc yields the posterior
% probabilities. Note that D*W*maxc yields the same classification results as
% the NN classifier.
%
% par = 'loo' - (default) compute leave-one-out optimization for a and b.
%               it is assumed that the first objects in the training set
%               constitute the representation set.
%       'all' - include all dissimilarities for optimization of k
%               (representation set should not be included in training set)
%
% See also mappings, datasets, classc, classd, testd, knnd
%
% $Id: discc.m,v 1.2 2001/07/10 16:35:58 pavel Exp $

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function W = discc(d,par)
			% set default leave-one-out
if nargin < 2
	par = 'loo';
end
			% empty call, to handle d*knnd, or d*knnd([],par)
if nargin < 1 | isempty(d)
	if isstr(par)	% untrained classifier, par = 'loo' or 'all'
		W = mapping('discc',par);
	else			% fixed classifier, par = knn
		W = mapping('discc','fixed',par);
	end
	return
end

[nlab,lablist,m,k,c,p,featlist] = dataset(d);
[clab,classlist] = renumlab(featlist);
c = size(classlist,1);
[cl,nc] = renumlab(classlist,lablist);
if size(nc,1) > c
	error('Object labels do not match representation set')
end

if isstr(par)
	
	% training (find parameters a and b)
			
	if strcmp(par,'loo')
		n = k*(k-1) + (m-k)*k;
		s = zeros(n,1);
		t = zeros(n,1);
		r = zeros(n,1);
		K = 0;
		for j=1:k
			if j == 1
				J = 2:k;
			elseif j < k
				J = [1:(j-1) j+1:k];
			else
				J = [1:k-1];
			end
			s(K+1:K+k-1) = d(j,J)';
			t(K+1:K+k-1) = repmat(nlab(j),k-1,1);
			r(K+1:K+k-1) = nlab(J);
			K = K+k-1;
		end
		for j=k+1:m
			s(K+1:K+k) = d(j,:)';
			t(K+1:K+k) = repmat(nlab(j),k,1);
			r(K+1:K+k) = nlab(1:k);
			K = K+k;
		end
		w = +(dataset(s,2-(t==r))*loglc);
		W = mapping('discc',w,lablist,k,c,1);
	elseif strcmp(par,'all')
		n = m*k;
		s = zeros(n,1);
		t = zeros(n,1);
		r = zeros(n,1);
		K = 0;	
		for j=1:m
			s(K+1:K+k) = d(j,:)';
			t(K+1:K+k) = repmat(nlab(j),k,1);
			r(K+1:K+k) = nlab(1:k);
			K = K+k;
		end
		w = +(dataset(s,2-(t==r))*loglc);
		W = mapping('discc',w,lablist,k,c,1);
	else
		error(['Unknown option ''' par ''''])
	end
	
else

	% testing for given mapping
		
	if isa(par,'mapping')
		w = +par;
		W = d*w(1) + w(2);
	else
		error('Second parameter should be string or mapping')
	end
end
