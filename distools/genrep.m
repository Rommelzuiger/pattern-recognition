%GENREP Generate a representation set
%
%   DR = GENREP(D,M,OPT)
%         OR
%   [DR,DT] = GENREP(D,M)
%
% INPUT
%   D   NxN dissimilarity dataset
%   M   Cx1 vector or scalar
%   OPT 'include' or 'exclude'
%
% OUTPUT
%   DR  - NxK dissimilarity dataset, if OPT = 'include' (DEFAULT);
%         K=CM or K=SUM(M)
%       - [N-K]xK dissimilarity dataset, if OPT = 'exlude';
%         K=CM or K=SUM(M)
%       - KxK dissimilarity dataset, if there is no third parameter;
%         K=CM or K=SUM(M)
%   DT  [N-K]xK remaining dissimilarity dataset, if there is no third
%       parameter; K=CM or K=SUM(M)
%
% DESCRIPTION
% If M is a scalar then a representation set of M objects per class is
% generated for the square (dis)similarity matrix D. D should be a dataset.
% DR has C*M columns if D has C classes. If M is a vector of length C,
% its elements determine the number of representation objects of the
% corresponding classes.
%
% If there is no third parameter, K=CM or K=SUM(M) objects are chosen
% for the represenation set. The resulting dissimilarity dataset DR is
% KxK and the remaining DT is [N-K]xK.
% If desired (OPT = 'exclude'), the objects used for the representation
% set are excluded from D. Per default, they are included.
%
% DEFAULT
% M   = 1
% OPT = 'include'
%

% Copyright: Robert Duin, r.duin@ieee.org, and
% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function [DR,DT] = genrep(D,N,str)

if nargin < 3,
  str = 'include';
end
if nargin < 2,
  N = 1;
end

nlab    = getnlab(D);
[m,k,c] = getsize(D);
labels  = getlab(D);
if size(N) == 1
  N = N*ones(1,c);
end
if length(N) ~= c
  error('Number of classes does not match the vector length.')
end

R = [];
for j = 1:c
  J = find(nlab==j);
  if length(J) < N(j)
    error('More objects requested than supplied.')
  end
  L = randperm(length(J));
  R = [R; J(L(1:N(j)))];
end

strl = lower(str);
if ~strcmp(strl,'exclude') & ~strcmp(strl,'include')
  error(['Unknown option demanded: ' str])
end

if nargout == 1
  L = [1:m]';
  if strcmp(strl,'exclude')
    L(R) = [];
  else
    L = [R; setdiff(L,R)];
  end
  DR = dataset(D(L,R),labels(L,:),'featlab',labels(R,:));
else
  if strcmp(strl,'exclude')
    error('EXCLUDE does not match here.')
  elseif strcmp(strl,'include')
    L = R;
    T = setdiff((1:m)',R);
  else
    ;
  end
  DR = dataset(D(L,R),labels(L,:),'featlab',labels(R,:));
  DT = dataset(D(T,R),labels(T,:),'featlab',labels(R,:));
end
return
