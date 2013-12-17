%SORT Sort in ascending order. Dataset overload

% $Id: sort.m,v 1.3 2007/03/22 07:45:54 duin Exp $

function [s,I] = sort(a,dim)
				
		nodatafile(a);
		
if nargin == 1
	[d,I] = sort(a.data);
else
	[d,I] = sort(a.data,dim);
end
s = setdata(a,d);
return
