%PE_KNNC K-Nearest Neighbor Classifier for PE spaces
%
%   [W,K,E] = PE_KNNC(A,K)
%   [W,K,E] = PE_KNNC(A)
%
% INPUT
%   A  PE dataset
%   K  Number of the nearest neighbors (optional; default: K is 
%      optimized with respect to the leave-one-out error on A)
%
% OUTPUT
%   W  k-NN classifier 
%   K  Number of the nearest neighbors used
%   E  The leave-one-out error of the KNNC
%
% DESCRIPTION  
% Computation of the K-nearest neighbor classifier for the PE dataset A. 
%
% Warning: class prior probabilities in A are neglected.
%
% SEE ALSO
% MAPPINGS, DATASETS, KNNC

% R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,k,e] = pe_knnc(a,k)

  if nargin < 2, k = []; end

	if nargin == 0 | isempty(a)
		w = mapping(mfilename,'untrained',{k});
		w = setname(w,'PE K-NN Classifier');
		
	elseif ~ismapping(k) % training
		
		if ~ispe_dataset(a)
			[w,k] = knnc(a,k);
		else
			if isempty(k)             % optimize k in PE space
				d = pe_distm(a);        % find PE distances
				[v,k,e] = knndc(d,k);   % use dis mat routine for optimisation k
      elseif nargout > 2
        e = testkd(pe_distm(a),k,'loo');
			end
			w = mapping(mfilename,'trained',{a,k},getlablist(a),size(a,2),getsize(a,3));
		end
		
	else % execution, testset is in a, trained mapping is in k
		
		%retrieve data
		trainset = getdata(k,1);
		k = getdata(k,2);
		d = pe_distm(a,trainset);
    [e,w] = testkd(d,k);         % confidences in w
    
  end
		
return