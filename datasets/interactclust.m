function o2c = interactclust(x,linkage,showfusion)

% ------------------------------------------------------------------
% o2c = interactclust(x,linkage)
%
% INPUT:
%     linkage: hierarchical linkage type
%              'a' average; 
%              's' single; 
%              'c': complete
%           x: data set: N X M (N samples, M features)
%  showfusion: (0/1) if set, fusion level plot is displayed, 
%              default is off
%
% OUTPUT:
%         o2c: object to cluster mapping
% ------------------------------------------------------------------

if ((nargin > 1) & (isa(x,'dataset'))),
  x = +x;
end;

if nargin < 3, showfusion = 0; end

fignum = 1;
dmeasure = 'e';
show_limit = 100;
NC = 16;

D = Clust_Distance(x,dmeasure);
Dend = Clust_Dendrogram(D,linkage);

% plot fusion level
if showfusion,
   tmp = Dend{2};
   merge_level = tmp(:,3);
   figure(fignum+2),clf,plot(merge_level(end:-1:end-min(length(merge_level)-1,15)),'.-');
   xlabel('Number of clusters','FontName','Verdana','FontSize',12)
   ylabel('Fusion level','FontName','Verdana','FontSize',12)
   grid
   set(gca,'FontName','Verdana','FontSize',12)
end
[obj2pos,minY] = Clust_Drawdend(fignum,Dend,'linear',show_limit);
set(fignum,'renderer','zbuffer')

% select number of clusters by clicking on dendrogram
exit_flag = 0;
ax2 = size(x,2);
Ns = size(x,1);
col = [];
o2c = [1:1:Ns];
while ~exit_flag,
   figure(fignum)
   [Nc_new,col_new,o2c_new,exit_flag] = Clust_Clickdend(Dend,'n',obj2pos,minY);
   if ~exit_flag,
      Nc = Nc_new;
      col = col_new;
      o2c = o2c_new;
      % create c2o from o2c
      c2o = [];
      for j = 1:Nc,
         c2o{j} = find(o2c == j);
      end
      for i=1:Nc,
         % select the color, to correspond to the dendrogram shading
         the_col = col(i,:);         
         if size(x,2) >= 2,
            if i==1,
               figure(fignum+1),clf;
            else
               figure(fignum+1),hold on;
            end
            sx = x(c2o{i},:);
            dots = plot(sx(:,1),sx(:,2),'.');
            set(dots,'color',the_col)
            if i==1,
               c = get(gcf,'children');
               set(c(1),'color','k')
            end
         end
      end
      set(gca,'FontName','Verdana','FontSize',12)   
   end
end
