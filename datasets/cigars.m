%CIGARS 200 objects with 2 features in 2 classes
%
%       X = CIGARS;
%
% The cigars dataset is an artificial dataset containing two very
% elongated but nonoverlapping classes. To find these two classes using
% an unsupervised clustering procedure is hard, and only suitable
% clustering methods will find the two clusters.
%
% SEE ALSO 
% DATASETS

function x = cigars

prdatasets(mfilename,1,'http://prtools.org/prcoursedata/cigars.mat');

user.desc = 'Artificial 2D dataset to test clustering procedures.';
user.link = '';

prload('cigars.mat');
x = setname(a,'Cigars');
x = setuser(x,user);

return
