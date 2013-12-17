%MTIMES Dataset overload of *
%
%    C = A*B
%
% This routines handles the dataset multiplication in case A or B is
% a dataset and none of these is a mapping. A or B may be a cell array
% of datasets as well as a scalar. 

% $Id: mtimes.m,v 1.4 2007/03/22 07:45:54 duin Exp $

function c = mtimes(a,b)
		
if isa(a,'double')
  if isscalar(a) 
    c = b*a;
  else % some matrix multiplication desired
    c = a*b.data;
  end
  return
elseif isa(a,'cell')
	if min(size(a)) ~= 1
		error('Only one-dimensional cell arrays are supported')
	end
	c = cell(size(a));
	for j=1:length(a)
		c{j} = a{j}*b;
	end
	return
end
if isa(b,'double')
  d = a.data*b;
  c = setdata(a,d);
  return
elseif isa(b,'cell')
	if min(size(b)) ~= 1
		error('Only one-dimensional cell arrays are supported')
	end
	c = cell(size(b));
	for j=1:length(b)
		c{j} = a*b{j};
	end
	return
end

d = a.data * b.data;
c = setdata(a,d,b.featlab);
return
