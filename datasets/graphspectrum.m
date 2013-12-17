%GRAPHSPECTRUM, Compute spectrum of a graph
%
%		C = GRAPHSPECTRUM(L,E)
%
% INPUT
%   L   Nx2 array with indices of all connected nodes in the graph, or
%       Nx3 array with indices of connected nodes and distances
%   E   Vector with N corresponding distances
%       Default: all distances equal to 1
%
% OUTPUT
%   S   Ranked set of eigenvalues
%
% DESCRIPTION
% The distance matrix over the graph is computed and used for 
% pseudo-Euclidean embedding. The related eigenvalues are ranked
% (largest values first) and returned.
%
% SEE ALSO
% GRAPH, KMST, GRAPHPATH, PLOTGRAPH, GRAPHDIST, PSEM

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands


function V = graphspectrum(L,ee);

	[n,p] = size(L);
	if p == 3
		e = L(:,3);
		L = L(:,[1 2]);
	elseif p == 2
		e = ones(n,1);
	else
		error('Input has wrong size')
	end
	
	if nargin > 1
		e = ee;
	end
	
	n = size(L,1);
	if nargin < 2
		e = ones(n,1);
	end
	d = graphdist(L);
	J = find(d==Inf);
	d(J) = zeros(size(J));
	[w,sig,V] = psem(d);
	V = V/sum(abs(V));
	if length(V)+1 < size(d,1)
		V = [V; zeros(size(d,1)-length(V)-1,1)];
	end
	V = fliplr(sort(V(:)'));

return
