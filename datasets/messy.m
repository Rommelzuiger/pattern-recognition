%MESSY 300 objects with 2 features in 2 clusters
%
%       X = MESSY;
%
% The messy dataset is an artificial dataset containing a very noisy
% two-cluster structure.
% Ideally, an clustering algorithm should be able to find two clusters here.
%
% SEE ALSO 
% DATASETS

function x = messy

prdatasets(mfilename,1,'http://prtools.org/prcoursedata/messy.mat');

user.desc = 'Artificial 2D dataset to test clustering procedures.';
user.link = '';

prload('messy.mat');
x = setname(a,'Messy');
x = setuser(x,user);

return
