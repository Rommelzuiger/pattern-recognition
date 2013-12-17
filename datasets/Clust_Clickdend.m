function [Nc,col,o2c,exit_flag] = Clust_Clickdend(Xd,scale,obj2pos,minY)

%---------------------------------------------------------------------------
% [Nc,col,o2c,exit_flag] = Clust_Clickdend(Xd,scale,obj2pos,minY)
%
% Select number of desired clusters by clicking on 
% the dendrogram
%
% INPUT:
%          Conn: connection history (see hc for details)
%         CHist: cluster history
%         scale: 'l' logscale on y-axis; 'n' normal scale on y-axis
%       obj2pos: obj2pos(k) is on-screen position of object k
%          minY: minimum Y-value on dendrogram plot
%
% OUTPUT:
%            Nc: number of clusters
%           col: colors used for dendrogram clusters
%           o2c: mapping from objects to clusters
%     exit_flag: set if right mouse button is used to select clustering
%                so that controlling program can exit
%---------------------------------------------------------------------------

Conn = Xd{2};
CHist = Xd{3};
Nobjects = size(CHist,1);

% pastel colors
new_map = make_own_colormap;
colormap('jet')
% get position of the cursor
x = ginput(1);
% left or right mouse button
butt = get(gcf,'selectiontype');
switch(butt),
case 'normal'
   exit_flag = 0;
case 'alt'
   exit_flag = 1;
   Nc = [];
   col = [];
   o2c = [];
   return
otherwise
   exit_flag = 0;
end
   
% determine the level from Conn
if x(2) < 0,
   x(2)  = eps;
end

Nc_at_cut = Nobjects - max(find(Conn(:,3)<x(2)));

% if clicked below is dendrogram bridge
if (isempty(Nc_at_cut)),
   % remove existing rectangles
   delete(findobj('tag','dend_rect'))
   % as many clusters as objects
   Nc_at_cut = Nobjects;
end
% get the sample cluster labels
[o2c,c2o] = Clust_Getobjects(Nc_at_cut, Xd);
Nc = Nc_at_cut;   
% construct the colormap based on the number of clusters
delta_color = floor(size(new_map,1)/Nc);
if (delta_color<1), delta_color=1; end
% select cluster colors uniformly spaced in palet
col_indeces = [1:delta_color:size(new_map,1)];
col = new_map(col_indeces,:);
% if previous rectangles exist, clear them
delete(findobj('tag','dend_rect'))
% draw rectangles including all elements in a cluster
% for all clusters at this level
for j=1:Nc,
   % elements in first cluster
   clust = find(o2c==j);
   % positions of these objects
   opos = obj2pos(clust);
   % object with smallest position
   [first_pos,ind] = min(opos);
   % object with largest position
   [last_pos,ind] = max(opos);
   % draw the rectangle
   R(j)=rectangle('Position',...
      [first_pos-1/5, minY,...
         last_pos-first_pos+2/5,...
         x(2)],'tag','dend_rect');
   % select the color
   the_col = col(j,:);
   % set face color of rectangle
   set(R(j),'facecolor',the_col);
   % set border color
   set(R(j),'edgecolor',the_col)
end
% Order the clusters according to size: largest to smallest
% Cluster 1 is largest, cluster Nobjects is smallest
% get sizes of the clusters to only plot 
% the largest clusters
L=[]; for i=1:Nc, L=[L;[i,length(c2o{i})]]; end
% sort on negative so that largest values ends up at top
L=sortrows(-L,2); L=-L;
% sort colors to correspond to the cluster order
col = col(L(:,1),:);
% new o2c
o2c = [];
for i = 1:Nc,
   obj_ind = c2o{L(i,1)}; 
   o2c(obj_ind) = i*ones(length(obj_ind),1);
end
























