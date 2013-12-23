%SETSIG Set signature for pseudo-Euclidean dataset or mapping
%
%   B = SETSIG(A,SIG)
%   B = A*SETSIG([],SIG);
%
% INPUT
%   A    Dataset or data matrix (doubles)
%   SIG  Signature, 2-element vector with numbers of positive
%        and negative dimensions, default: all positive
%
% OUTPUT
%   B   PE dataset
%
% DESCRIPTION
% Distances in a PE dataset are computed in a special way, see PE_DISTM
% using the signature of the PE space. This routine stores the signature
% in the user-field of the dataset, so it can be used by any routine that
% needs it. The signature is determined during a pseudo-Euclidean mapping,
% e.g. by PSEM. Datasets applied to such a mapping obtain the signature
% automatically. The signature is a vector of two elements equal to the
% numbers of positive and negative eigenvalues found during the process of
% PE mapping. Their sum should thereby be equal to the dimensionality of
% the space, i.e. the number of features in A.
%
% Note that V = W*SETSIG converts a pseudo-Euclidean mapping W into a
% mapping V to the associated space.
%
% SEE ALSO
% DATASETS, MAPPINGS, PE_EM, PE_DISTM, GETSIG

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = setsig(a,sig)

  if nargin < 2, sig = []; end
  if nargin < 1 | isempty(a)
    b = mapping(mfilename,'fixed',{sig});
    b = setname(b,'Set signature');
    return
  end
  
  if isempty(sig), sig = [size(a,2) 0]; end
  if length(sig) ~= 2
    error('The PE signature should consist of a vector of two elements');
  end
  a = dataset(a);
  if sum(sig) ~= size(a,2)
    error('The signature does not match the feature size')
  end
  b = setuser(a,sig,'pe_signature');
  
return
  
  
  