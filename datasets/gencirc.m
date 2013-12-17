%
% GENCIRC  Generate circular data.
% 
%  DATA = GENCIRC (N) generates N samples drawn from a Gaussian distribution
%  around a circle. Default: N = 100.
 
function data = gencirc (n)

	if (nargin < 1), n = 100; end;

	radii = ones(n,1) + 0.1 * randn(n,1);
	angle = (0:(2*pi/n):2*pi*(n-1)/n)';

	data(:,1) = radii .* sin(angle);
	data(:,2) = radii .* cos(angle);

	labs = ones(n,1);

	data = prdataset(data,labs);

return
