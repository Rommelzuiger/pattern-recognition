%GETLABLISTNAMES Get the names of all label lists
%
%	NAMES = GETLABLISTNAMES(A)
%
% INPUT
%   A      - Dataset
%
% OUTPUT
%   NAMES  - Character array with label list names
%
% DESCRIPTION
% All label list names are returned in NAMES. In case no multiple
% labeling is set, NAMES = 'default'.
%
% SEE ALSO
% DATASETS, MULTI_LABELING, ADDLABELS, ADDLABLIST, CHANGELABLIST,
% CURLABLIST, SETLABLISTNAMES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function names = getlablistnames(a)
				
	if ~iscell(a.lablist)
		names = 'default';
	else
		names = a.lablist{end,1};
	end
	
return
