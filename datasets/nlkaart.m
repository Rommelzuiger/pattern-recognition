% NLKAART  Show a map of the Netherlands
%
%  NLKAART (F) shows the map of the Netherlands in figure F 
%  (default: current figure).

function nlkaart (f)

	if (nargin < 1), f = gcf; end; 

	figure (gcf);
	clf; hold on;

	[im,cmap] = imread ('steden.bmp');
	%im = flipud (255-9*double(im));
	im = flipud(im);
	imshow(im,cmap);

	axis image off;

return;
