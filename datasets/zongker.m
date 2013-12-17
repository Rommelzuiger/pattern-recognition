%ZONGKER Distances between digits by Zongker
%
%       X = ZONGKER;
%
% The Zongker dataset is a distance matrix between digits, where the
% distance between images of the digits are measured by warping the
% images of the distances, and by measuring how much effort it takes.
%
% SEE ALSO 
% DATASETS

function x = hall

prdatasets(mfilename,1,'http://prtools.org/prcoursedata/dist_digit_zongker.mat');

user.desc = 'Distance matrix between digits by Zongker.';
user.link = '';

prload('dist_digit_zongker.mat');
x = setname(a,'Zongker');
x = setuser(x,user);

return
