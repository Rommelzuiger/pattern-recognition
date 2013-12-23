%PE_PARZENC Parzen Classifier for PE spaces
%
%   [W,H,E] = PE_PARZENC(A,H)
%   [W,H,E] = PE_PARZENC(A)
%
% INPUT
%   A  PE dataset
%   H  Smoothing parameter (optional; default: H is optimized
%      with respect to the leave-one-out error on A)
%
% OUTPUT
%   W  Parzen classifier 
%   H  Number of the nearest neighbors used
%   E  The leave-one-out error
%
% DESCRIPTION  
% Computation of the Parzen classifier for the PE dataset A. 
%
% Warning: class prior probabilities in A are neglected.
%
% SEE ALSO
% MAPPINGS, DATASETS, PARZENC, PARZENDDC

% R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,h,e] = pe_parzenc(a,h)

  if nargin < 2, h = []; end

  if nargin == 0 | isempty(a)
    w = mapping(mfilename,'untrained',{h});
    w = setname(w,'PE Parzen Classifier');
    
  elseif ~ismapping(h) % training
    
    if ~ispedataset(a)
      [w,h] = parzenc(a,h);
    else
      if isempty(h)             % optimize h for PE space
        d = sqrt(pe_distm(a));  % find PE distances
        [v,h] = parzenddc(d,h); % use dis mat routine for optimisation h
      end
      if nargout > 2
        e = testpd(sqrt(pe_distm(a)),h,'loo');
      end
      w = mapping(mfilename,'trained',{a,h},getlablist(a),size(a,2),getsize(a,3));
    end
    
  else % execution, testset is in a, trained mapping is in h
    
    %retrieve data
    trainset = getdata(h,1);
    h = getdata(h,2);
    
    d = sqrt(pe_distm(a,trainset));
    [e,w] = testpd(d,h);         % confidences in w
    
  end
    
return