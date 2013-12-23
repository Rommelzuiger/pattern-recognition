%TESTPD Test Parzen classifier for dissimilarity data
%
%   [E,C] = TESTPD(D,H,PAR)
%
% INPUT
%   D    Dissimilarity dataset. Object labels are assumed to be the
%        true labels. Feature labels are assumed to be the labels of
%        the objects they are related to.
%   H    Desired smoothing parameter.
%   PAR  'loo' - leave-one-out option. This should be used if
%        the objects are related to themselves. If D is not square,
%        it is assumed that the first sets of objects in columns and
%        rows match.
%        'all' - use all objects (default).
%
% OUTPUT
%   E    Estimated error
%   C    Matrix with confidences, size M x N, if D has size M x L and
%        the labels are given for N classes. 
%
% SEE ALSO
% DATASETS, KNND, TESTKD

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [e,F] = testkd(d,h,par)

if nargin < 3 | isempty(par), par = 'all'; end
if nargin < 2 | isempty(h), [w,h] = parzenddc(d); end

isdataset(d);

nlab     = getnlab(d);
lablist  = getlablist(d);
featlist = getfeat(d);
[m,k,c]  = getsize(d);
p        = getprior(d);

[clab,classlist] = renumlab(featlist);
c = max(clab);
%[cl,nc,labl] = renumlab(classlist,lablist);
%if size(labl,1) > c
%	error('Object labels do not match representation set')
%end
                   % correct for different classlist - lablist orders
J = matchlablist(classlist,lablist);
classlist = lablist;
clab = J(clab);

if strcmp(par,'loo')
	% get rid of leave-one-out problems
	km = min(k,m);
	dmax=max(max(+d))*100;
	d(1:km,1:km) = d(1:km,1:km) + dmax*eye(km);
elseif ~strcmp(par,'all')
	error(['Unknown option ''' par ''''])
end

s = exp(-(d.^2)/(2*h*h));

F = zeros(m,c);
for j=1:c
  L = find(clab==j);
  F(:,j) = mean(s(:,L),2)*p(j);
  if strcmp(par,'loo')
    if length(L) == 1
      F(L,j) = 0;
    else
      F(L,j) = length(L)*F(L,j)/(length(L)-1); % correct for LOO
    end
  end
end
%F = F ./ (sum(F,2)*ones(1,c));
F = setdata(d,F,classlist);
e = F*testc;