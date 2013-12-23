%GETSIG Retrieve the signature from a pseudo-Euclidean dataset or mapping
%
%   SIG = GETSIG(W)
%   SIG = GETSIG(A)
%
% INPUT
%   W    PE mapping, W = PE_EM(D), if D is a dissimilarity matrix
%   A    Dataset, vectors in PE space, A = D*W
%
% OUTPUT 
%   SIG  Signature, 2-element vector with numbers of 
%        positive and negative dimensions
%
% SEE ALSO
% DATASETS, MAPPINGS, SETSIG, PE_EM, PE_DISTM

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function sig = getsig(a)

  if isdataset(a)
    sig = getuser(a,'pe_signature');
    if isempty(sig)
      sig = [size(a,2) 0];
    end
  elseif ismapping(a)
    if ispsem(a)
      a = pe_em(a);
    end
    if ispe_em(a)
      sig = getdata(a,'sig');
    else
      sig = [size(a,2) 0];
    end
  else % doubles
    sig = [size(a,2) 0];
  end
    
return
  
  
  