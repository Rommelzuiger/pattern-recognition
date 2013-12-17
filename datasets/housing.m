%HOUSING 506 objects with 13 features in two classes
%
%       X = HOUSING;
%
% The housing dataset is real-world dataset contains 13 fields of
% information on areas in the suburbs of Boston. The goal is to predict
% whether the median price of a house in each area is larger than or
% smaller than $20,000. 
%
% SEE ALSO 
% DATASETS

function x = housing

prdatasets(mfilename,1,'http://prtools.org/prcoursedata/prhousing.mat');

user.desc = 'The housing dataset from UCI. The goal is to predict whether the median price of a house in each area is larger than or smaller than $20,000. ';
user.link = 'http://www.ics.uci.edu/~mlearn/MLRepository.html';
fl = { 'crime    ';
'large    ';
'industry ';
'river    ';
'nox      ';
'rooms    ';
'age      ';
'work     ';
'highway  ';
'tax      ';
'education';
'aa       ';
'status   '};
prload('prhousing.mat');
x = setname(a,'housing');
x = setfeatlab(x,fl);
x = setuser(x,user);

return
