%GETLABLIST Get label list of dataset
%
%   LABLIST = GETLABLIST(A,STRING)
%
% Returns the label list of a dataset A.
% If STRING equals 'string' the label list is converted to characters,
% which may be useful for display purposes and annotation of graphs.
%
% SEE ALSO MULTI_LABELING

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: getlablist.m,v 1.4 2006/09/26 13:07:34 duin Exp $

function lablist = getlablist(a,string)
	
	if iscell(a.lablist)
		n = a.lablist{end,2};
		lablist = a.lablist{n,1};
	else
		lablist = a.lablist;
	end
	
	if (nargin > 1) & (string == 1 | strcmp(string,'string')) & (~isstr(lablist))
		lablist = [repmat('Class ',length(lablist),1) num2str(lablist)];
	end
	
return;
