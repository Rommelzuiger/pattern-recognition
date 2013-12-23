%DMSTSPM Find the Shortest Paths along K Minimum Spanning Trees
%
%   W = DMSTSPM (D,K,P)
%
% INPUT
%   D   NxN Symmetric dissimilarity matrix or dataset
%   K   Number of minimum spanning trees; (optional; default: 1)
%   P   Object number, P = 1..N;
%       (optional; default: P with the largest distances to all other objects)
%
% OUTPUT
%   W   Mapping that determines the shortest paths
%
% DESCRIPTION
% Determines the shortest paths along a neighborhood graph defined by K minimum
% spanning trees. MSTs are found by the Prim's algorithm, starting from the
% point P. The result, however, does not depend on P. New vertices are projected
% on an edge of the existing neighborhood graph such that the distances
% between a new object and the vertices of this edge are the smallest. Then,
% the distances along the shortest paths are computed.
%
% DEFAULT
% K = 1
% P = max(dist)
%
% SEE ALSO
% KMST, DSPATH, DSPATHS

% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty of EWI, Delft University of Technology, The Netherlands and 
% School of Computer Science, University of Manchester




function W = dmstspm(D,K,p)

if nargin < 3 | isempty(p),
  [ss,p]= max(sum(+D));
end
if nargin < 2 | isempty(K),
  K = 1;
end
if nargin < 1 | isempty(D)
   W = mapping(mfilename,{K,p});
   W = setname(W,'dmstspm');
   return
end


if isdataset(D) | isa(D,'double'),

  if ismapping(K),
    n_orig = getsize_in(K);
    pars   = getdata(K);

    L   = pars{1};    % list of edges in the MSTs
    dm  = pars{2};    % weights of the edges
    K   = pars{3};    % number of MSTs

    W   = D;
    D   = +D;
    [m,n] = size(D);
    if n ~= n_orig,
      error('The number of columns in the test data does not match.')
    end

    is_porig = (m > n &  sum(diag(+D(1:n,:)))<=1e-16);
    is_orig  = (m == n & sum(diag(+D))<=1e-16);

    LL  = [L; L(:,[2 1])];
    ddm = [dm; dm];
    DD  = inf * ones(n_orig);
    DD((LL(:,2)-1)*n_orig+LL(:,1)) = ddm;
    DD(1:n_orig+1:end)=0;

    if is_orig | is_porig,        % if (a part of) test data is the original distance matrix
      Dp = dspaths(DD);
    end

    if ~is_orig,
      if is_porig,
        pos = n;
        Dp  = [Dp; zeros(m-n,n)];
      else
        pos = 0;
        Dp  = zeros(m,n);
      end

      dmn = zeros(m-pos,2);
      dd  = zeros(m-pos,1);
      Ln  = zeros(m-pos,2);
      D = D(pos+1:m,:);

      [dmin,I] = min(D,[],2);
      LL = L(:);
      for i=1:length(I)
        J = find(LL == I(i));
        [dd(i),z] = min(ddm(J));
        Ln(i,1:2) = [I(i) LL(J(z))];
      end

      % Use the cosine law to determine the distances of the points
      % projected onto the closest edge
      X   = (dmin.^2 + dd.^2 - D(Ln(I,2)*(m-pos)+I).^2) ./(2*dd);
      dmn = [X dd-X];
      Z   = find(X < 0);  % make corrections if distances are negative
      if ~isempty(Z)
        dmn(Z,:) = [-X(Z) dd(Z)-X(Z)];
      end
      Z = find(X > dd); % make corrections if distances are too large
      if ~isempty(Z)
        dmn(Z,:) = [X(Z) -dd(Z)+X(Z)];
      end

      for i=1:length(I)
        dnew = inf*ones(1,n);
        dnew(Ln(i,:)) = dmn(i,:);
        Q = dspaths([DD dnew'; dnew 0]);  % Floyd's algorithm is used, since it is fast in Matlab
        Dp(pos+i,:) = Q(n+1,1:n);
      end
    end


    if isdataset(W),
      W = setdata(W,Dp);
    else
      W = Dp;
    end
    return;

  else

    D = +D;
    [n,m] = size(D);
    if n ~= m | ~issym(D,1e-12),
      error('D should be symmetric.')
    end
    D = 0.5*(D+D');

    [Lmst,d] = kdmst(D,K,p);
    W = mapping(mfilename,'trained',{Lmst,d,K},[],n,n);
    W = setname(W,'dmstspm');
  end
end
return;
