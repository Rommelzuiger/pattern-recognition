% BOOMERANGS  Generation of boomerang-shaped classes
%   A = BOOMERANGS(N) generates N samples from two boomerang-shaped
%   classes. N can be a vector with the number of samples per class.

function a = boomerangs (n)

	if (nargin < 1), error ('Specify number of samples/class'); end;

	p = [0.5 0.5];
	n = genclass(n,p);

	t = pi * (-0.5 + rand(n(1),1));

	xa = 0.5*cos(t)           + 0.025*randn(n(1),1);
	ya = 0.5*sin(t)           + 0.025*randn(n(1),1);
	za = sin(2*xa).*cos(2*ya) + 0.025*randn(n(1),1);

	t = pi * (0.5 + rand(n(2),1));

	xb = 0.25 + 0.5*cos(t)    + 0.025*randn(n(2),1);
	yb = 0.50 + 0.5*sin(t)    + 0.025*randn(n(2),1);
	zb = sin(2*xb).*cos(2*yb) + 0.025*randn(n(2),1);

	a = prdataset ([xa ya za; xb yb zb], genlab(n));

return
