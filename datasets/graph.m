%GRAPH  Construct graph from coordinates by Delaunay tessellation
%
%	  L = GRAPH(X)
%
% INPUT
%   X     MxK matrix with coordinates of N points in K dimensions
%
% OUTPUT
%   L     Nx3 array with indices of connected nodes and their distances
%         (edge length)
%
% SEE ALSO
% PLOTGRAPH, DISTGRAPH, GRAPHPATH, KMST

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function L = graph(x,varargin)

	T = delaunayn(x);
	[m,k] = size(T);
	p = nchoosek(size(T,2),2);
	L = zeros(m*p,2);
	n = 1;
	for j=1:m
		S = combnk(T(j,:),2);
		L(n:n+size(S,1)-1,:) = S;
		n = n+p;
	end
	L = unique(L,'rows');
	d = (x(L(:,1),:)-x(L(:,2),:)).^2;
	d = sqrt(sum(d,2));
	L = [L,d];
