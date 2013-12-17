function rbwhite_map = rbwhite

% Generates the yellow and blue colormap with no expression set to zero
% number of points per up-ramp of a color
Np_up = 49;

x_up = [0:1:Np_up];

y_up = x_up/max(x_up);
U0 = y_up;
U1 = y_up(1:end-1);
U2 = y_up(2:end);

y_down = -sort(-x_up)/max(x_up);
D0 = y_down;
D1 = y_down(2:end);
D2 = y_down(1:end-1);

Red = [ones(1,Np_up),y_down];
Green = [y_up(1:end-1),y_down];
Blue = [y_up, ones(1,Np_up)];

%rbgreen_map=[R',W',G'];  % odd map

rbwhite_map = [Red',Green',Blue'];
