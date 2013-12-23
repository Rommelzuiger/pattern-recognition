%ISSQUARE Test on square dissimilarity matrix
%
%   OK = ISSQUARE(D)
%   ISSQUARE(D)
%
% INPUT
%   D      Dataset
%
% OUTPUT
%   OK     1 if the matrix D is square and 0, otherwise.
%
% DESCRIPTION
% True is D is a square dissimilarity matrix dataset. This includes
% the check whether feature labels equal object labels.
% If called without an output argument ISSQUARE generates an error
% if D is not square.

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function OK = issquare(d)
isdataset(d);
[m,k] = size(d);

if m == k
  n  = nlabcmp(getfeatlab(d),getlabels(d));
  OK = (n == 0);
else
  OK = 0;
end

if nargout == 0 & OK == 0
  error([newline '---- Square dissimilarity matrix expected ----'])
end
