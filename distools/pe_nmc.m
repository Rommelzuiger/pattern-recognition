%PE_NMC Nearest Mean Classifier for PE spaces
%
%   W = PE_NMC(A)
%
% INPUT
%   A  PE dataset
%
% OUTPUT
%   W  Nearest mean classifier 
%
% DESCRIPTION  
% Computation of the nearest mran classifier for the PE dataset A. 
%
% Warning: class prior probabilities in A are neglected.
%
% SEE ALSO
% MAPPINGS, DATASETS, NMC

% R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = pe_nmc(a,w)

  if nargin == 0 | isempty(a)
    w = mapping(mfilename,'untrained');
    w = setname(w,'PE Nearest Mean');
    
  else
    
    if ~ispe_dataset(a)
      w = nmc(a);
    else
      u = meancov(a);
      sig = getsig(a);
      u = setsig(u,sig);
      w = pe_knnc(u,1);
    end
    
  end
    
return