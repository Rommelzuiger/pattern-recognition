%DISCHECK Dissimilarity Matrix Check
%
%   N = DISCHECK(D)
%         OR
%   N = DISCHECK(D,LABLIST,FLAG)
%
% INPUT
%   D       NxN dissimilarity matrix or dataset
%   LABLIST List of labels
%   FLAG    0 - check square dissimilarity matrix
%           1 - check square similarity matrix
%
% OUTPUT
%   N      1/0 if D is/isnot a proper dissimilarity matrix
%
% DESCRIPTION
% Returns error if 
% 	- D is not a square dissimilarity matrix with feature labels 
% 		equal to object labels or feature labels of D are not in LABLIST.
% 	- FLAG==0 & there are non-zero elements on the diagonal of D
% 	- FLAG==0 & there are negative values in D
%

% Copyright: R.P.W. Duin, r.duin@ieee.org
% and Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester

function out = discheck(D,lablist,flag)

if nargin < 3, flag = 0; end
if nargin < 2 | isempty(lablist)
  [m,k,c]  = getsize(D);
  lablist  = getlablist(D);
  featlist = getfeatlab(D);
  if m ~= k
    if nargout == 0
      error('Square (dis)similarity matrix D is expected.')
    else
      out = 0;
      return
    end
  end
else
  featlist = getfeat(D);
  c = size(lablist,1);
end
[nn,nf,fl] = renumlab(lablist,featlist);
if max(nf) > c
  if nargout == 0
    error('Feature labels do not match the class labels.')
  else
    out = 0;
    return
  end
end

if any(+D < 0) & ~flag,
  if nargout == 0
    error('Dissimilarities should be nonnegative.')
  else
    out = 0;
    return
  end
end

dd = diag(+D);
if any(dd ~= 0) & ~flag,
  if max(dd) < 1e-5 & min(dd) >= 0
     D(1:m+1:end) = 0;
     prwarning(3,'The diagonal has small non-zero values. They are set to zero.'); 
  else
    if nargout == 0
      error('The diagonal of D should be a zero vector.')
    else
      out = 0;
      return
    end
  end
end

if nargout > 0
  out = 1;
end

return;
