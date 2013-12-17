%GRAPHPATH, Compute shortest paths in a graph
%
%		[PATH,D] = GRAPHPATH(N,L,E)
%
% INPUT
%   N   Index of start node
%   L   Nx2 array with indices of all connected nodes in the graph, or
%       Nx3 array with indices of connected nodes and distances
%   E   Vector with N corresponding distances
%       Default: all distances equal to 1
%
% OUTPUT
%   PATH Cell array with paths from node N to all other nodes
%   D    Cell array with edge length of paths in PATH
%
% DESCRIPTION
% The shortest paths are found from node N to all other nodes 
% in the graph defined by {L,E}
%
% SEE ALSO
% GRAPH, GRAPHDIST, KMST

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [nodepath,pathdist] = graphpath(nstart,L,ee);

[n,p] = size(L);
if p == 3
	e = L(:,3);
	L = L(:,[1 2]);
elseif p == 2
	e = ones(n,1);
else
	error('Input has wrong size')
end

if nargin > 2
	e = ee;
end

L = [L; fliplr(L)];
e = [e(:);e(:)];
k = max(L(:));

nodedist = repmat(inf,1,k);  % distance of node to source
nodeconn = cell(1,k);        % connections per node
nodepath = cell(1,k);        % path from source to each node
pathdist = cell(1,k);        % distances along path
edgelen  = cell(1,k);
for j=1:k
	J = find(L(:,1)==j);
	nodeconn{j} = L(J,2);
	edgelen{j} = e(J);
end

K = nstart;                  % active nodes
nodepath{nstart} = nstart;
nodedist(nstart) = 0;
pathdist{nstart} = 0;

while ~isempty(K)
	Knew = [];
	for j=1:length(K)
		C = nodeconn{K(j)};         % connecting nodes
		for i=1:length(C)
			ndist = nodedist(K(j)) + edgelen{K(j)}(i);
			if ndist < nodedist(C(i))
				nodedist(C(i)) = ndist;
				nodepath{C(i)} = [nodepath{K(j)} C(i)];
				pathdist{C(i)} = [pathdist{K(j)} edgelen{K(j)}(i)];
				Knew = [Knew C(i)];
			end
		end
	end
	K = Knew;
end