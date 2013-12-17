%HALL 140 objects with 2 features
%
%       X = HALL;
%
% The hall dataset is an artificial dataset with a very clear clustering
% structure. You may discuss if the dataset contains three, or 14
% clusters....
%
% SEE ALSO 
% DATASETS

function x = hall

prdatasets(mfilename,1,'http://prtools.org/prcoursedata/hall.mat');

user.desc = 'Artificial 2D dataset to test clustering procedures.';
user.link = '';

prload('hall.mat');
x = setname(a,'Hall');
x = setuser(x,user);

return
