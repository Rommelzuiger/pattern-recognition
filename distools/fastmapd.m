%FASTMAPD Fast Linear Map of Euclidean distances
%
% 	W = FASTMAPD(D,K)
%
% INPUT
%   D   NxN symmetric dissimilarity matrix (dataset)
%   K   Chosen dimensionality (optional; K is automatically detected)
%
% OUTPUT
%   W   Linear map of Euclidean distances D into a Euclidean space
%
% DEFAULT
%   K   = []
%
% DESCRIPTION
% This is an implementation of the FastMap algorithm. D is assumed to be 
% a square Euclidean distance matrix and K is the desired dimensionality of 
% the projection. However, if the dimensionality M of the original data is
% smaller than K, then only a mapping to a M-dimensional space is returned.
%
% If D does not have a Euclidean behavior, a projection is done to the fewest
% dimensions that explain the Euclidean part the best, i.e. do not deal with
% negative squared distances.
%
% LITERATURE
% 1. C. Faloutsos and K.-I. Lin, FastMap: A Fast Algorithm for Indexing, 
% Data-Mining and Visualization of Traditional and Multimedia Datasets, 
% ACM SIGMOD, International Conference on Management of Data, 163-174, 1995.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [W,P] = fastmapd(D,K)

if nargin < 2, 
	K = []; 
end
if nargin < 1 | isempty(D)
   W = mapping(mfilename,K); 
   W = setname(W,'Fastmapd');
   return
end


if (isdataset(D) | isa(D,'double')) 
  if ismapping(K),
    % APPLY THE MAPPING
    pars = getdata(K);
    P    = pars{1};
    dd   = pars{2}; 
    X    = pars{3}; 
    K    = pars{4}; 

		nlab = getnlab(D);
		n    = size(D,1);
    D    = +D.^2;
    if all(diag(D) == 0) & issym(D) & length(X) == n,
      Y = X;
    else
      for k=1:K
        Y(:,k) = (D(:,P(k,1)) + dd(k)*ones(n,1) - D(:,P(k,2)))./(2*sqrt(dd(k)));
        D      = (D - distm(Y(:,k),X(:,k)));
      end
		end		
    W = dataset(Y,nlab);
    return
  end
end





% TRAIN THE MAPPING
if isdataset(D) 
	nlab = getnlab(D);
end
discheck(D);

n = size(D,1);
D = +D.^2;
if isempty(K), 
	K = n-1; 
end

k      = 0;     % Current dimensionality
Z      = 1:n; 	% Index of points
tol    = 1e-12; % Tolerance when the distances become close to zero 
tolNE  = 1e-10; % Tolerance for the non-Euclidean behavior
NEucl  = 0; 		% Non-Euclidean behavior
finito = 0;

while ~finito,
  k = k+1;
  [mm,O] = max((D(Z,Z)),[],2);
  [mm,i] = max(mm);
 	P(k,1) = Z(O(i));		% P is the index of pivots
  P(k,2) = Z(i);
  Z([O(i) i]) = [];
	dd(k)  = D(P(k,1),P(k,2));

  if dd(k) >= tol,
    X(:,k) = (D(:,P(k,1)) + dd(k)*ones(n,1) - D(:,P(k,2)))./(2*sqrt(dd(k)));
    D      = (D - distm(X(:,k)));    
	else
		P(k,:) = [];
	end

	if ~NEucl,
		NEucl = NEucl | any(D(:) < -tolNE);
		if NEucl,
			KKe = k;
			KK  = k;
		end 
	end

  finito = (k >= K) | (dd(k) <= tol) | NEucl; 
end
if ~NEucl,
	KK = k;
end

if NEucl,
  disp('Distances do not have a Euclidean behavior!');
	disp(['The Euclidean part of the distances is explained in ' num2str(KKe) 'D.']);
else
	if KK < K, 
  	disp(['The distances are perfectly explained in ' num2str(KK) 'D.']);
	end
end

W = mapping(mfilename,'trained',{P,dd,X,KK},[],n,KK);
W = setname(W,'Fastmapd');
return



