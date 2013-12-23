%NHGRAPH Find a neighborhood graph and its shortest paths
%
%   [NG,L,DL,D] = NHGRAPH(A,OPT,PAR)
%
% INPUT
%   A   NxN Weight (dissimilarity) matrix/dataset representing a weighted graph
%   OPT Option method (optional; default: 'MST')
%       'NN'     - based on nearest neighbor neighborhoods
%       'MST'    - based on minimum spanning trees (MSTs)
%       'EPS'    - based on eps-neighborhoods
%       'MST-NN' - based on a single MST followed by nearest neighbor neighborhoods
%       'MST-EPS'- based on a single MST followed by eps-neighborhoods
%   PAR Parameter (optional; default: 3 for 'NN' and 'MST-NN'; 1 for 'MST';
%       0.1*avr(W) for 'EPS' and 'MST-EPS')
%       Integer K, 1<=K<N, for 'NN', 'MST-NN'
%       Integer K, 1<=K<N, for 'MST'
%       Real positive value for 'EPS', 'MST-EPS'
%
% OUTPUT
%   NG  NxN Weight matrix (dataset) representing the neighborhood graph
%   L   List of edges in the neighborhood graph
%   DL  List of edge weights of the neighborhood graph
%   D   NxN matrix of shortest paths
%
% DESCRIPTION
% Finds a neighborhood graph for the graph described by the weight matrix A.
% A is NxN matrix representing a graph G(V,E). V is a set of vertices, |V| = N,
% and E is a set of edges. If there is no edge between i and j then A(i,j) = INF.
%
% DEFAULT
%   M   = 'MST'
%   PAR = 3 ('NN','MST-NN'), 1 ('MST') or 0.1*avr(W) ('EPS','MST-EPS')
%
% SEE ALSO
% KMST, DSPATH, DSPATHS

% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty of EWI, Delft University of Technology, The Netherlands and
% School of Computer Science, University of Manchester 


function [V,L,DL,D] = nhgraph(A,opt,par)

if nargin < 2 | isempty(opt),
  opt = 'MST';
end

opt = upper(opt);
if nargin < 3 | isempty(par),
  switch opt
    case 'MST', par = 1;
    case {'MST','NN','MST-NN'}, par = 3;
    case {'EPS','MST-EPS'}, par = 0.1*mean(A(:));
    otherwise
			error(['Wrong method', opt]);
  end
end


V = +A;
[n,m] = size(V);
if n~= m,
  error('Weight matrix A should be square.');
end


L  = [];
DL = [];
switch opt
  case 'MST'
    [L,DL,V] = kmst(V,par);

  case {'NN', 'MST-NN'}
    K = par;
    if K ~= round(K) | K < 2 | K >= n,
      error ('K should be the number of nearest neighbors.');
    end
    if strcmp(opt,'MST-NN'),
      [L,DL,V] = kmst(V,1);
    end
    [VV,I] = sort(V);
    for i=1:n
      V(i,I((2+K):end,i)) = inf;
    end
    if nargout > 1,
      for i=1:n
        L  = [L; i*ones(K,1) I(2:1+K,i)];
        DL = [DL; VV(i,I(2:1+K,i))'];
      end
    end

  case {'EPS','MST-EPS'}
    if strcmp(opt,'MST-EPS')
      [L,DL,V] = kmst(V,1);
    end
    deps = par;
    if deps <= 0,
      error ('EPS should be positive.');
    end
    warning off;
    V = V./(V <= deps);
    V = min(V,inf);
    if nargout > 1,
      [I J] = find(V < inf);
      L  = [L; I J];
      DL = [DL; V(J*(n-1)+I)];
    end
    warning on;
  otherwise
    error('Unknown method.')
end

LL = [L; L(:,[2 1])];
V  = inf*ones(n,m);
V((LL(:,2)-1)*n+LL(:,1)) = +A((LL(:,2)-1)*n+LL(:,1));

if isdataset(A),
  V = setdata(A,V);
end

if nargout > 3,
  D = dspaths(V);
end
return;
