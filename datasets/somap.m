function somap(a,nunits,lattice,shape)

% SOMAP a routine based on the somtoolbox
% (http://www.cis.hut.fi/projects/somtoolbox/)
% and PRtools for visualizing datasets.
%
% SOMAP(A,NUNITS,LATTICE,SHAPE)
%
% It trains a SOM with a 2D input map on the 
% given dataset (A) employing NUNITS input 
% nodes and the specified LATTICE and SHAPE.
%
% Possible alues for parameters:
%   LATTICE: 'hexa' (default),'rect' 
%     SHAPE: 'sheet' (default),'cyl','toroid'}
% 
% It employs the somtoolbox visualization 
% tools to display the positions of the input 
% nodes in the dataspace as well as the U-matrix
% (distance between nodes in data space) and then
% also projects the SOM in the dataspace in a 
% gridded scatterd plot.

%---extract dataset feature labels from dataset 
%---and stick in SOM struct
lab = getfeatlab(a);
clab = [];
for i=1:size(lab,1),
    clab{i} = lab(i,:);
end
sData = som_data_struct(+a,...
    'name',getname(a),...
    'comp_names',clab);

sD = som_normalize(sData,'var');
if nargin < 2,
    sMap = som_make(sD);
elseif nargin < 3,
    sMap = som_make(sD, 'munits', nunits);
elseif nargin < 4,
    sMap = som_make(sD, 'munits', nunits, lattice);
else
    sMap = som_make(sD, 'munits', nunits, lattice, shape);
end

% sMap = som_make(sD, 'munits', nunits, lattice, shape);
sMap = som_autolabel(sMap,sD,'vote');

%---do the standard SOM display
figure(1),clf
colormap('rbwhite')
som_show(sMap,'norm','d')

if 0,
    %---Calculate the histograms: number of obejects of each class
    %---in a given cell
    labs = getnlab(a);
    lablist = unique(labs);
    h = [];
    for i=1:size(lablist,1),
        h = [h,som_hits(sMap,sD.data(labs==i,:))];
    end
    cols = [0 0 0;
        0 1 0;
        1 1 0;
        1 0 1;
        0 1 1;
        1 0 1;
        0 0 1;
        1 0 0;
        1 1 1];
    som_show_add('hit',h,'MarkerColor',cols(1:size(h,2),:),'Subplot',1)
end

%---display the grids in scatterd - gridded format
figure(2),clf
M = som_denormalize(sMap.codebook,sMap);
D = som_denormalize(sD.data,sD); 
scatterd(a,'gridded')
k=0;
[n,p] = size(getdata(a));
for i=1:p,
    for j=1:p,
        k=k+1;
        subplot(p,p,k);
        hold on;
        som_grid(sMap,...
            'Coord',M(:,[j,i]),...
            'LineWidth',1,...
            'LineColor','k',...
            'MarkerSize',2);
    end;
end
return


