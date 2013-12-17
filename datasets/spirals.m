%SPIRALS 194 objects with 2 features in two classes
%
%       X = SPIRALS;
%
% The spirals dataset is an artificial dataset with two classes that are
% intertwined in a spirally way:-)
%
% SEE ALSO 
% DATASETS

function x = spirals

prdatasets(mfilename,1,'http://prtools.org/prcoursedata/spirals.mat');

user.desc = 'Artificial 2D dataset to test classification procedures.';
user.link = '';

prload('spirals.mat');
x = setname(a,'Spirals');
x = setuser(x,user);

return
