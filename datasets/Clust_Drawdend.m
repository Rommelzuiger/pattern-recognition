function [obj2pos,minY] = Clust_Drawdend(figure_index,dend,scale,show_limit)

% ------------------------------------------------------------------------------------
% [obj2pos,minY] = Clust_Drawdend(figure_index,dend,scale)
%
% INPUT:
%  figure_index: matlab figure index where dendrogram should be rendered
%          Conn: connection history (see hc for details)
%      Nobjects: Number of objects on data set
%         scale: 'log' logscale on y-axis; 'linear' normal scale on y-axis
%    show_limit: if more objects in dendrogram than show_limit: do not render labels
%
% OUTPUT:
%       obj2pos: obj2pos(k) is on-screen position of object k
%          minY: minimum Y-value on dendrogram plot
% ------------------------------------------------------------------------------------

lincol = 'k';

Cl = dend{1};
Conn = dend{2};
CHist = dend{3};
Nobjects = size(CHist,2);

figure(figure_index),clf

% number of levels in dendrogram
K = Nobjects - 1;

% logarithmic scale on Y-axis
switch scale, 
 case 'log', 
  % logscale on Y-axis 
  set(gca,'YScale','log')
  % we need a y-position reference for the log scale
  % since log(0) = -inf
  minY = min(Conn(:,3));
 case 'linear', 
  minY = 0;
  set(gca,'YScale','linear');
 otherwise 
  error('Unknown scale')
end

% position objects equidistantly along the x-axis 
% at position y = 0 in the order preserved in Cl
object_x_pos = [1:1:Nobjects];

% mapping of objects to position
for i = 1:Nobjects, obj2pos(Cl(i)) = object_x_pos(i); end

% remove the standard ticklabels
set(gca,'xticklabel',[])
% remove tickmarks
set(gca,'xtick',[])

% (x, y) coordinates where vertical dendrogram line should start
join_pos_X = obj2pos;
join_pos_Y = minY * ones(1, Nobjects);

% actual plot loop
% (TOTAL LENGTH OF LONGEST DENDRAGRAM LEG IS CLUSTER DISTANCE)
for i = 1:K,

 % get elements to join eL and eR
 eL = Conn(i, 1);
 eR = Conn(i, 2);

 % max of the join_pos_Y of the two being joined
 highest_join_pos = max(join_pos_Y(eL), join_pos_Y(eR));

 % new y_pos
 new_join_pos_Y = Conn(i, 3);

 % vertical line for left-most object
 line([join_pos_X(eL), join_pos_X(eL)],[join_pos_Y(eL), new_join_pos_Y],[.5 .5],...
       'Color',lincol);

 % vertical line for right-most object
 line([join_pos_X(eR), join_pos_X(eR)],[join_pos_Y(eR), new_join_pos_Y],[.5 .5],...
       'Color',lincol);

 % update the join_pos_Y of the left most element
 join_pos_Y(eL) = new_join_pos_Y;

 % draw the "bridge"
 line([join_pos_X(eL), join_pos_X(eR)],[join_pos_Y(eL), join_pos_Y(eL)],[.5 .5],...
       'Color',lincol);

 % update the join_pos_X of the left most element
 join_pos_X(eL) =  (join_pos_X(eL) + join_pos_X(eR)) / 2;

end;

% scale axis to fit it all ...
Ax = axis;
axis([0, max(object_x_pos) + 1, minY, Ax(4)]);
set(gca,'FontName','Verdana','FontSize',12)

% only show the labels if fewer than 200
if Nobjects < show_limit,
   switch scale 
   case 'log',
      % determine position of text
      del = (log10(Ax(4))-log10(minY))/60;  	% 60 equal increments in linear scale
      tpos = 10^(log10(minY)-del);	       	% translate to log scale: one del below min 
      % labels of objects
      for i = 1:Nobjects, 
         text(object_x_pos(i), tpos, sprintf('%d',double(Cl(i))),...
            'HorizontalAlignment','left',...
            'VerticalAlignment','middle',...
            'Rotation',270,...
            'FontName','Verdana',...
            'FontWeight','Normal',...
            'FontSize',10);
      end;
   case 'linear',			% labels of objects
      for i = 1:Nobjects, 
         text(object_x_pos(i), -Ax(4)/60, sprintf('%d',double(Cl(i))),...
            'HorizontalAlignment','left',...
            'VerticalAlignment','middle',...
            'Rotation',270,...
            'FontName','Verdana',...
            'FontWeight','Normal',...
            'FontSize',10);
      end;
   otherwise,
 		error(['EUSError: Unknown Scale ' scale '.']);
    end;
end






