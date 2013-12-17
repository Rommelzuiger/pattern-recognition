% GENDATDD Generate difficult classes
%
%  A = GENDATDD (N,D) generates a 2-class data set containing N points in 
%  a D-dimensional distribution. In this distribution, 2 randomly chosen 
%  dimensions contain fully separated Gaussian distributed data; the others 
%  contain unit covariance Gaussian noise.
%
%  Defaults: N = 100, D = 2.

function data = gendatdd (n, d)

	if (nargin < 2), d = 10;  end;
	if (nargin < 1), n = 100; end;

	data = +gauss(n,2*ones(1,d),5*eye(d));

	a = +gauss(floor(n/2),[0 0],[3 -2.5; -2.5 3]);
	b = +gauss(ceil(n/2), [4 4],[3 -2.5; -2.5 3]);

	p = randperm(d);
	data(:,p(1:2)) = [a; b];

	labs = [ ones(floor(n/2),1); 2*ones(ceil(n/2),1) ];
	
	data = prdataset(+data,labs);

return

