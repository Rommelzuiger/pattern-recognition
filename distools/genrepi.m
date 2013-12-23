%GENREPI Generate indices for representation, learning and testing sets
%
% 	[R,L,T] = GENREPI(D,K,M)
%
% INPUT
%   D   NxN dissimilarity dataset
%   K   Number of objects for the representation set R
%   M   Number of objects for the learning set L; M >= K
%
% OUTPUT
%   R  Indices for the representation set
%   L  Indices for the learning set
%   T  Indices for the test set
%
% DESCRIPTION
% Random selection of M indices for a representation set R and N indices for 
% a learning set L. L contains R (M <= N). Indices of the remaining objects 
% are stored in T.
%
% DEFAULT
% M   = N
%

% Copyright: Robert Duin, r.duin@ieee.org, and
% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester
%

function [R,L,T] = genrepi(D,k,m)
if nargin < 3, 
	m = k; 
end

[nr,nc,c] = getsize(D);
if nr ~= nc,
	error('Dissimilarity matrix should be square.');	
end

if m > nr | k > nr
	error('K and M should not exceed the total number of objects in D.')
end

if k > m
	error('K should not be larger than M.');
end
		
R = []; 
L = []; 
T = [];
C = classsizes(D);

for j = 1:c
	J  = findnlab(D,j); 
	LL = randperm(C(j)); 
	LL = LL(1:m); 
	RR = LL(1:k);
	R  = [R;J(RR)]; 
	L  = [L;J(LL)];
end
T = setdiff([1:nr]',L);
	
