%TESTKD Test k-NN classifier for dissimilarity data
%
%   [E,C] = TESTKD(D,K,PAR)
%
% INPUT
%   D    Dissimilarity dataset. Object labels are assumed to be the
%        true labels. Feature labels are assumed to be the labels of
%        the objects they are related to.
%   K    Desired number of neighbors to take into account, default K = 1.
%   PAR  'LOO' - leave-one-out option. This should be used if
%          the objects are related to themselves. If D is not square,
%          it is assumed that the first sets of objects in columns and
%          rows match.
%        'ALL' - use all objects (default).
%
% OUTPUT
%   E    Estimated error
%   C    Matrix with confidences, size M x N, if D has size M x L and
%        the labels are given for N classes. Note that for K < 3 these
%        confidences are derived from the nearest neigbor distances
%        and that for K >= 3 they are the Bayes estimators of the 
%        neighborhood class probabilities.
%
% DESCRIPTION
% TESTK is based on just counting errors and does not weight with 
% testobject priors.
%
% SEE ALSO
% DATASETS, KNNDC

% Copyright: R.P.W. Duin, r.duin@ieee.org
% Faculty EWI, Delft University of Technology


function [e,F] = testkd(d,knn,par)

if nargin < 3 | isempty(par), par = 'all'; end
if nargin < 2 | isempty(knn), knn = 1; end

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
%	error('Object labels do not match representation set.')
%end
                   % correct for different classlist - lablist orders
J = matchlablist(classlist,lablist);
classlist = lablist;
clab = J(clab);

if strcmp(upper(par),'LOO')
	% get rid of leave-one-out problems
	km = min(k,m);
	dmax=max(max(+d))*2;
	d(1:km,1:km) = d(1:km,1:km) + dmax*eye(km);
elseif ~strcmp(upper(par),'ALL')
	error(['Unknown option ''' par ''''])
end


% find class frequencies in representation set			
% r = zeros(1,c);
% for j=1:c
% 	r(j) = length(find(clab==j));
% end

%D = ones(m,c);
[DD,L] = sort(+d',1);			% sort distances
L = clab(L);
for j = 1:c     				% find label frequencies
	F(:,j) = sum(L(1:knn,:)==j,1)';
end
K = max(F');
%for j = 1:c
%	K = min(K,r(j));  
%	J = reshape(find(L==j),r(j),m); % find the distances to the
%	J = J(K+[0:m-1]*r(j));		% objects of that neighbor
%	D(:,j) = DD(J)';			% number for all classes
%end
                                % estimate posterior probabilities
%if knn > 2					    % use Bayes estimators on frequencies
	F = (F+1)/(knn+c);
	%else						% use distances
%  F = sigm(log(sum(D,2)*ones(1,c)./(D+realmin) - 1 + realmin));
%end
F = F ./ (sum(F,2)*ones(1,c));
F = setdata(d,F,classlist);
%e = F*testc; % goes wrong in case of LOO testing (empty classes)
labf = F*labeld;  
e = nlabcmp(labf,getlabels(F))/m;