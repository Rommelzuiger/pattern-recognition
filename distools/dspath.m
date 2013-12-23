%DSPATH  Single Shortest Path in a (Dissimilarity) Graph
%
%   [D,P] = DSPATH (A,i,j)
%
% INPUT
%   A   NxN Weight / dissimilarity matrix or dataset
%   i,j Vertices for which the shortest path should be found
%
% OUTPUT
%   D   Distance of the shortest path
%   P   List of edges of the shortest path
%
% DESCRIPTION
% Determines the shortest path between two vertices i and j. The Dijkstra's
% algorithm is used. Currently, it is not optimized fro speed in Matlab.
% A is a NxN matrix of weights (or a distance matrix) representing a graph
% G(V,E). V is a set of vertices, |V| = N, and E is a set of edges. If there
% is no edge between i and j then A(i,j) = INF. Use DSPATHS(A,'F') if all
% shortest paths should be computed.
%
% SEE ALSO
% DSPATHS
%
% LITERATURE
% Any book on graph algorithms or basic algorithms in computer science.

% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty of EWI, Delft University of Technology, The Netherlands and 
% School of Computer Science, University of Manchester


function [d, pt] = dspath(A,s,t)

A = +A;
[n,m] = size(A);

I = (1:n);
dmin(I)  = 1e10;
final(I) = 0;
pred(I)  = -1;


I(s)     = [];
dmin(s)  = 0;
final(s) = 1;
last     = s;

while ~final(t)
  dminnew = dmin(last) + A(last,I);
  Z = find(dminnew < dmin(I));
  dmin(I(Z)) = dminnew(Z);
  pred(I(Z)) = last;
  [ss,last]  = min(dmin(I));
  last = I(last);
  final(last) = 1;
  I = setdiff(I,last);
end

if nargout > 1,
  if pred(t) ~= -1   % there is a path
    k  = t;
    pt = t;
    while k ~= s,
      pt = [pred(k) pt];
      k  = pred(k);
    end
  else
    pt = [];
  end
end

d = dmin(t);
return;
