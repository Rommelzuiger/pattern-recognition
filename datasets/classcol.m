%CLASSCOL Display class colors
%
%	CLASSCOL(LABLIST)
%	CLASSCOL(W)
%
% Displays an image with the colors as used by PLOTC(W) in a SCATTERD plot.

function classcol(arg)
if nargin < 1
	error('No argument given')
end
if ismapping(arg)
	lablist = getlabels(arg);
else
	lablist = arg;
end

c = size(lablist,1)
lablist
figure
image(repmat([1:c]',1,1));
colormap(0.5+hsv(10)*0.5);
set(gca,'xtick',[]);
set(gca,'ytick',[1:c]);
set(gca,'yticklabel',lablist);

