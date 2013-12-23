%DRSSCC Dissimilarity-based Random Subspace Combining Classifier 
%
% 	W = DRSSCC(D,V,N,M,CRULE)
%
% INPUT
%   D 		NxK Dissimilarity dataset
% 	V 		Base untrained classifier
% 	N 		Desired number of representation objects per class
% 	M 		Number of base classifiers to be generated
% 	CRULE Combining rule
%
% OUTPUT
%   W 		Trained classifier
%
% DESCRIPTION
%
% DEFAULT
% 	V 		= NMSC
% 	N 		= 2
% 	M 		= 11
% 	CRULE = VOTEC   
%

% Copyright: R.P.W. Duin, duin@ieee.org and
% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function W = drsscc(D,v,N,M,crule)

if nargin < 5, 
	crule = votec; 
end
if nargin < 4, 
	M = 11; 
end
if nargin < 3, 
	N = 2; 
end
if nargin < 2, 
	v = nmsc; 
end
if nargin < 1 | isempty(D)
  W = mapping('DRSSCC',{v,N,M,crule});
  return
end


nlab   = getnlab(D);
[m,k]  = size(D);
[flab,flist] = renumlab(getfeat(D));
c = size(flist,1);

if length(N) == 1,
	N = N*ones(c,1);
elseif length(N) ~= c,	
	error('N should be either a scalar or a vector with the length C.');
else
	;	
end	

W = [];
for j = 1:M
  R = genreps(flab,N);
  W = [W D*(cmapm(k,R)*v)];
end
W = traincc(D,W,crule);



	
function R = genreps(nlab,N)
% Generate n objects per class
c = max(nlab);
R = [];
for j = 1:c
  J = find(nlab==j);
  L = randperm(length(J))';
  R = [R; J(L(1:N(j)))];
end
