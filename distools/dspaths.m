%DSPATHS  All Shortest Paths in a (Dissimilarity) Graph
%
%   D = DSPATHS (A,ALG)
%
% INPUT
%   A   NxN Weight/dissimilarity matrix or dataset
%   ALG Algorithm used (optional; default: 'F')
%       'F' Floyd's algorithm; fast in Matlab
%       'D' Dijkstra's algorithm; slow in Matlab
%
% OUTPUT
%   D   NxN matrix of the shortest paths
%
% DESCRIPTION
% Determines all shortest paths in a (dissimilarity) graph. A is a NxN matrix
% of weights (or a dissimilarity matrix) representing a graph G(V,E). V is
% a set of vertices, |V| = N, and E is a set of edges. If there is no edge
% between i and j then A(i,j) = INF. D is the matrix of the values of shortest
% paths between all pairs of vertices. Either Floyd's or Dijkstra's algorithm
% is used.
%
% DEFAULT
% ALG = 'F'
%
% SEE ALSO
% DSPATH
%
% LITERATURE
% Any book on graph algorithms or basic algorithms in computer science.

% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty of EWI, Delft University of Technology, The Netherlands and 
% School of Computer Science, University of Manchester

function A = dspaths (A,alg)

if nargin < 2,
  alg = 'F';
else
  alg = 'D';
end

V = +A;
[n,m] = size(V);
if n~=m,
  error('Weight matrix W should be square.')
end

switch upper(alg)
  case 'F',
    for k=1:n
      V = min(V,repmat(V(k,:),n,1)+repmat(V(:,k),1,n));
    end
  case 'D',
    VV = zeros(n,n);
    for i=1:n
      for j=1:n
        VV(i,j) = dspath(V,i,j);
      end
    end
    V = VV;
otherwise
  error('Wrong algorithm.')
end

if isdataset(A),
  A = setdata(A,V);
else
  A = V;
end
