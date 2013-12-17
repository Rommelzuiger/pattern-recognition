function [x,y] = genregres(nrx,noise)
%GENREGRES Generate regression data
%
%     Y = GENREGRES(NRX,SIGMA)
%
% Generate an artificial regression dataset [X,Y] with:
%
%    y = sin(4x) + noise. 
%
% where noise is Gaussian distributed with standard deviation sigma.
%
%  Y = GENREGRES(100)
% generates 100 x,y pairs with data and default noise (sigma = 0.1).
%
%  x = (0:0.01:1)';
%  y = genregres(x,0);
% generates the true function along the x-axis, with zero noise.
%
if nargin<2
  noise = 0.1;
end

if (length(nrx)>1)
  x = nrx;
  nrx = size(x);
  y = sin(4*x) + noise*randn(nrx);
else
  x = rand(nrx,1);
  y = sin(4*x) + noise*randn(nrx,1);
end

if nargout<2
	x = gendatr(x,y);
end

return

