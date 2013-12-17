function my_map = make_own_colormap

% number of points per up-ramp of a color
Np_up = 50;

x_up = [0:1:Np_up];
y_up = x_up/max(x_up);
y_down = (max(x_up(2:Np_up))-x_up(2:Np_up))/max(x_up(2:Np_up));
y = [y_up,y_down];

Z=zeros(1,Np_up);
O=ones(1,Np_up);
U = y_up;
D = y_down;

B = [O,D,Z,U];
G = [U,O,D,Z];
R = [Z,U,O,D];

my_map=0.5.*[R',G',B']+0.5;

if 1>2,
colormap(my_map)

figure(2), clf
colormap(my_map)
rgbplot(my_map)
colorbar

% plays with colormap to enhance cluster analysis
figure(1),clf
colormap(my_map);
color_index = 1;
%new_map=map(17:61,:);
new_map=my_map;
delta_color = floor(length(new_map)/Nc);
for i = 1:Nc,
 subplot(round(sqrt(Nc)),ceil(sqrt(Nc)),i)
 set(gca,'color',new_map(color_index,:));
 line([0,1],[0,1],'color','k')
 color_index = color_index + delta_color;
 title(color_index)
end
end