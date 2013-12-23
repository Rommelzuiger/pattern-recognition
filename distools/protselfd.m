%PROTSELFD Forward Prototype Selection for Dissimilarity Matrices
%
%   [W,E,KOPT] = PROTSELFD(D,K,PAR)
%    W  = D*PROTSELFD([],K,PAR)
%
% INPUT
%   D     Dataset, square dissimilarity matrix
%   K     Integer, desired number of prototypes
%   PAR  'LOO' - leave-one-out option. This should be used if
%          the objects are related to themselves. If D is not square,
%          it is assumed that the first sets of objects in columns and
%          rows match.
%        'ALL' - use all objects (default).
%
% OUTPUT
%   W     Selection mapping ('feature selection')
%   E     Error stimate as a function of number of selected prototypes
%         (only reliable for prototype sizes >= class size)
%   KOPT  Estimate for best size in avoiding peaking
%
% DESCRIPTION
% This procedure for optimizing the representation set of a 
% dissimilarity matrix is based on a greedy, forward selection of
% prototypes using the leave-one-out error estimate of the 1NN rule
% as a criterion. As this is computed on the given distances in
% D, the procedure is based on sorting and counting only and is
% thereby fast.
%
% Note that the search continues untill K prototypes are found.
% This might be larger than desired due to peaking (curse of
% dimensionality, overtraining). Therefor an estimate for the
% optimal number of prototype is returned in KOPT. 
%
% The prototype selection may be applied by C = B*W(:,1:KSEL),
% in which B is a dissimilarity matrix based on the same 
% representation set as A (e.g. A itself) and C is a resulting 
% dissimilarity matrix in which the KSEL (e.g. KOPT) best prototypes
% are selected.
%
% REFERENCE
% E. Pekalska, R.P.W. Duin, and P. Paclik, Prototype selection for
% dissimilarity-based classification, Pattern Recognition, 
% vol. 39, no. 2, 2006, 189-208.
%
% SEE ALSO
% KNNDC, DISEX_PROTSELFD

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% 

function [R,e,D,J,nlab,clab] = protselfd(D,ksel,par,J,e,nlab,clab)

if nargin < 2, ksel = []; end
if nargin < 3 | isempty(par), par = 'all'; end

if nargin < 4 % user call
  
  if nargin < 1 | isempty(D)  % allow for D*protselfd([],pars)
    R = mapping(mfilename,'untrained',{ksel,par});
    R = setname(R,'Forward Prototype Sel');
    return
  end
  
  [m,k,c] = getsize(D);
  if isempty(ksel), ksel = k; end
  if strcmp(par,'loo') | strcmp(par,'LOO')
    if k > m
      error('More rows than columns expected for dissimilarity matrix')
    end
    discheck(D(1:k,:));
    D(1:k,:) = D(1:k,:) + 1e100*eye(k); % get rid of diagonal for LOO
  end
  
  %Initialise by the centre of the largest class
  cc = classsizes(D);
  [cmax,n] = max(cc); % n is the largest class
  lablist = getlablist(D);
  nlab = getnlab(D);
  clab = renumlab(getfeatlab(D),lablist);
  R = find(nlab == n);
  C = find(clab == n);
  dd = +D(R,C);
  [dmin,rmin] = sort(dd,1); % find one but most remote object
  [dmin,cmin] = min(dmin(end-1,:)); % find prototype for which this is minimum
  R = C(cmin);
  
  e = zeros(1,ksel);
  [nlab,clab] = renumlab(getlabels(D),getfeatlab(D));
  [dd,J] = min(+D(:,R),[],2);
  e(1) = sum(clab(R(J)) ~= nlab);
  
  % this will be a deep recursive call !!!
  prwaitbar(ksel,'Forward prototype selection')
  [R,e,D,J,nlab,clab] = protselfd(D,ksel,R,J,e,nlab,clab);
  prwaitbar(0);
  e = e(1:length(+R))/m;
  R = featsel(k,R);
  
  % Find optimal number of prototypes in avoiding peaking
  
  Jopt = find(e==min(e));
  D = floor((Jopt(end)+Jopt(1))/2);
  
  % done!
  
else  % internal call
  
  R = par;
  [m,k,c] = getsize(D);
  d = +D;
  S = [1:k];
  S(R) = [];
  emin = inf;
  dmin = inf;
  r = length(R);
  prwaitbar(ksel,r);
  for j=S
    [ds,n] = min([d(m*(R(J')'-1)+[1:m]'),d(:,j)],[],2);
    cclab = [clab(R(J)') repmat(clab(j),m,1)];
    ee = sum(cclab(m*(n-1)+[1:m]') ~= nlab);
    de = sum(ds);
    if ee < emin | ((ee == emin) & (de < dmin))
      emin = ee;
      jmin = j;
      JJ = [J repmat(r+1,m,1)];
      Jmin = JJ(m*(n-1)+[1:m]');
      Rmin = [R jmin];
      dmin = de;
    end
  end

  if emin <= e(r) | 1
    e(r+1) = emin;
    R = Rmin;
    if (r+1) < ksel
	[R,e,D,J,nlab,clab] = protselfd(D,ksel,R,Jmin,e,nlab,clab);
    end
  end
  
end
