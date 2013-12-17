%PLOTGRAPH Plot a graph
%
%   PLOTGRAPH(A,L,S)
%
% INPUT
%   A   Nx2 or Nx3 Data matrix
%   L   List of edges belonging to the graph
%
% DESCRIPTION
% Plots data and edges of a graph. Applicable only for 2D or 3D data.
% If A is unknown, but dissimilarity data are given, then a 2D or 3D 
% approximation can be found by PSEM.
%
% SEE ALSO
% KMST, PSEM, GRAPHDIST, GRAPHPATH

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function plotgraph (A,L,S)

if nargin < 3, S = 'k'; end
A = +A;
[n,m] = size(A);

len = length(L);
K   = round (len/(n-1));

if m == 2,
  plot(A(:,1),A(:,2),'.r')
  text(A(:,1),A(:,2),num2str([1:n]'));
  grid on;
  hold on;
  for i=1:len
    plot([A(L(i,1),1) A(L(i,2),1)],[A(L(i,1),2) A(L(i,2),2)],S);
  end

elseif m == 3,

  plot3(A(:,1),A(:,2),A(:,3),'.r')
  text(A(:,1),A(:,2),A(:,3),num2str([1:n]'));
  grid on;
  hold on;
  for i=1:len
    plot3([A(L(i,1),1) A(L(i,2),1)],[A(L(i,1),2) A(L(i,2),2)],[A(L(i,1),3) A(L(i,2),3)],S);
  end
else
  error('Plotting is only possible for two or three variables.');
  return;
end

hold off
if K == 1,
  title('Minimum Spanning Tree');
else
  title([num2str(K) ' Minimum Spanning Trees']);
end
return;
