%TRICLUST 140 objects with 2 features
%
%       X = TRICLUST;
%
% The triclust dataset is a simple artificial dataset with a very clear
% clustering structure with three clusters. Any sensible clustering
% algorithm should find these clusters...
%
% SEE ALSO 
% DATASETS

function x = triclust

prdatasets(mfilename,1,'http://prtools.org/prcoursedata/triclust.mat');

user.desc = 'Artificial 2D dataset to test clustering procedures.';
user.link = '';

prload('triclust.mat');
x = setname(a,'Triclust');
x = setuser(x,user);

return
