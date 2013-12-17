%PROSTATE Load prostate data
%
%    X = PROSTATE
%
function x = prostate

prdatasets(mfilename,1,'http://prtools.org/prcoursedata/prostate.mat');

prload prostate;
a = +a;
y = a(:,end);
x = gendatr(a(:,1:(end-1)),y);

fl = {'lcoval' 'lweight' 'age' 'lbph' 'svi' 'lcp' 'gleason' 'pgg45'};
x = setfeatlab(x,fl);

return
