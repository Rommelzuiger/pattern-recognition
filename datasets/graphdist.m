%GRAPHDIST, Compute distances in a graph
%
%		D = GRAPHDIST(L)
%
% INPUT
%   L   Nx2 array with indices of connected nodes, or
%       Nx3 array with indices of connected nodes and distances
%
% OUTPUT
%   D   Full square distance matrix
%
% DESCRIPTION
% The distances between all nodes in the graph are computed by
% adding all distances of the connecting nodes. Unconnected
% nodes have distance INF.
%
% SEE ALSO
% GRAPH, KMST, GRAPHPATH, PLOTGRAPH

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands


function g = graphdist(L);

[n,p] = size(L);
if p == 3
	e = L(:,3);
	L = L(:,[1 2]);
elseif p == 2
	e = ones(n,1);
else
	error('Input has wrong size')
end
R = 1e10/max(e);
e = round(e*R);   % reduce accuracy to avoid bit ripling
k = max(L(:));    % size of distance matrix (highest node index)
                  % construct initial distance matrix
g = repmat(inf,k,k);
g(1:k+1:k*k) = zeros(1,k);
                  % substitute given distances
for j=1:n
	g(L(j,1),L(j,2)) = e(j);
end
h = min(g,g');
                  % take care that edges are given two-way
L = [L; fliplr(L)];
x = [e(:);e(:)];
                  % loop as long as changes
while any(h(:) ~= g(:))
	g = h;
	for j=1:length(x)
		h(L(j,1),:) = min(h(L(j,1),:),h(L(j,2),:)+x(j));
	end
	h = min(h,h');
end

g = h/R;           % reset scale