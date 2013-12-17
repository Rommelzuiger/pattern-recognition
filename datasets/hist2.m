% HIST2  2D histogram.
%
%   Z = HIST2(Y) bins the elements of the data Y into 10x10 equally spaced
%   containers and returns the number of elements in each container. 
%
%   Z = HIST2(Y,M), where M is a scalar, uses M bins.
% 
%   [Z,X,Y] = HIST2(...) also returns the position of the bin centers in 
%   X and Y.
% 
%   HIST2(...) without output arguments produces a histogram bar plot of
%   the results.

function [z, x, y] = hist2 (data, n, linespec)

	if (nargin < 3), linespec = 'b';                 end;
	if (nargin < 2), n = 10;                         end;
  if (nargin < 1), error ('not enough arguments'); end;

	if isdataset(data)
		data = +data;
	end;

	if (size (data,2) ~= 2)
		error ('data should be an n x 2 data set');
	end;

	a = data(:,1); 
  b = data(:,2); 

	mn = min(min(min(a)),min(min(b))); 
	mx = max(max(max(a)),max(max(b))); 
	st = (mx-mn)/n;

	x_axis = mn+st/2:st:mx-st/2; y_axis = x_axis;
	[x,y] = meshgrid (x_axis, y_axis);

	for i = 1:n

		ia = find (a  >= mn+(i-1)*st); ta = a(ia);  tb = b(ia);
		ia = find (ta <  mn+i*st);     ta = ta(ia); tb = tb(ia);

		for j = 1:n

  		ib = find (tb  >= mn+(j-1)*st); ttb = tb(ib); 
  		ib = find (ttb <  mn+j*st);     ttb = ttb(ib);

			z(j,i) = length(ttb);

		end;

	end;

	if (nargout < 1)
		my_bar3 (y_axis, z, linespec); 
		set (gca, 'YDir', 'normal');
	end;

return
