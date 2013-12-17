% LINES5D  Generates three 5-dimensional lines
%
%   DATA = LINES5D(N);
%
%   Generates a data set of N points, on 3 non-crossing, non-parallel lines
%   in 5 dimensions. Default: N = 100.

function data = lines5d (n)

	n1 = floor(n/3);
	n2 = n-2*n1;

  s1 = [0 0 0 1 0];
  s2 = [1 1 1 0 0];
  s3 = [0 1 0 1 0];
  s4 = [1 1 1 1 1];
  s5 = [0 1 1 0 1];
  s6 = [1 0 1 1 1];
  c1 = [0:1/(n1-1):1]';
  c2 = [0:1/(n2-1):1]';
  a  = c1*s1 + (1-c1)*s2;
  a  = [a; c1*s3 + (1-c1)*s4];
  a  = [a; c2*s5 + (1-c2)*s6];

	data = prdataset(a,[ones(n1,1); 2*ones(n1,1); 3*ones(n2,1)]);

return
