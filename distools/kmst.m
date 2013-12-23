%KMST Finds K Minimum Spanning Trees based on a Distance Matrix
%
%   [L,DL,D] = KMST(D,K,P)
%
% INPUT
%   D   NxN Symmetric dissimilarity matrix  or dataset
%   K   Number of minimum spanning trees (optional, default: 1)
%   P   Object number, P = 1..N (optional, default: object P with the
%       largest distances to all other objects)
%
% OUTPUT
%   L   K(N-1)x2 List of edges in the minimum spanning tree
%   DL  K(N-1)x1 List of the corresponding minimal distances
%   D   NxN Dissimilarity matrix (dataset) with removed edges of L
%       (set to INF)
% DEFAULT
%   K = 1
%   P = Object number which has the largest distances to all other objects
%
% DESCRIPTION
% Finds K minimum spanning trees (MSTs) based on the given symmetric distance
% matrix D. D(i,j)=INF and D(j,i)=INF denote the lack of edge (connection)
% between i and j. Starting from the first MST, the procedure finds the
% current MST and then removes all the edges from D. It repeats until K MSTs
% are determined. The search for MSTs starts from the object P. The result does
% not depend on P, however, the ordering does.
%
% The implementation relies on the Prim's algorithm.
%
% SEE ALSO
% PLOTGRAPH, DISTGRAPH, GRAPHPATH

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function [L,d,D] = kmst(D,k,p)

if nargin < 3,
  [ss,p]= max(sum(+D));
end

if nargin < 2,
  k = 1;
end

[n,m] = size(D);
if n ~= m,
  error('D should be a square matrix.');
end
prec = 1e-12;
if ~issym(D,prec),
  prwarning(1,'D should be symmetric. It is made so by averaging.');
end
D = 0.5*(D+D');

if p < 0 | p > n |  p ~= round(p),
  error(['P should be an integer in [1,' int2str(n) ']']);
end
if k < 0 | k > n | k ~= round(k),
  error(['K should be an integer in [1,' int2str(n) ']']);
end

isdset = isdataset(D);

if isdset,
  DD = D;
  D  = +D;
end
L  = [];
d  = [];
for z=1:k
  Lmst = zeros(n-1,2);
  dmst = zeros(n-1,1);
  Z    = setdiff(1:n,p);
  [dmst(1),z] = min(D(p,Z));
  Lmst(1,:)   = [p Z(z)];
  for i=2:n-1
    LL   = Lmst(1:i-1,:);
    Jmst = unique(LL(:))';
    J    = setdiff(1:n,Jmst);
    [Dmin,Z]    = min(D(Jmst,J),[],2);
    [dmst(i),k] = min(Dmin);
    Lmst(i,:)   = [Jmst(k) J(Z(k))];
  end

  LL = [Lmst; Lmst(:,[2 1])];
  %for u=1:length(LL), D(LL(u,1),LL(u,2)) = inf; end
  D((LL(:,2)-1)*n+LL(:,1)) = inf;

  L = [L; Lmst];
  d = [d; dmst];
end

if isdset,
  D = setdata(DD,D);
end
return;
