% CIRCLES3D  Create a data set containing 2 circles in 3 dimensions.
%
%   DATA = CIRCLES3D(N) creates a data set containing N points (default: 100).

function data = circles3d (n)

	if (nargin < 1), n = 100; end;

	n2a = floor(n/2);
	n2b = ceil(n/2);
	ha = 0:(2*pi/n2a):2*pi*(n2a/(n2a+1)); ha = ha';
	hb = 0:(2*pi/n2b):2*pi*(n2b/(n2b+1)); hb = hb';

	a = [ sin(ha) cos(ha) zeros(n2a,1) ];
	b = [ sin(hb) cos(hb) ones(n2b,1)  ];

	data = prdataset ([a;b],[ones(n2a,1); 2*ones(n2b,1)]);

return
