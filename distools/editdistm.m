%EDITDISTM  Edit Distance Matrix between Strings
%
%   D = EDITDISTM (A,B,COST,SC)
%
% INPUT
%   A     Cell structure of N strings
%   B     Cell structure of M strings (optional; default: B = A)
%   COST  Vector of edit costs, COST = [cost_ins cost_del cost_sub]
%         (optional; default: all cost equal 1, COST = [1 1 1])
%   SC    Scale (1) or not (0) the computed distance to [0,1] by
%         max {K*COST(3)+(L-K)*COST(1), K*cost(2)+L*COST(1),...
%         L*COST(2)+K*COST(1), K*cost(3)+(L-K)*COST(2)]}, where
%         K =length(A{i}) and L = length(B{j}) and K > L.
%         (optional; default: 0, i.e. no scaling)
%
% OUTPUT
%   D     NxM Dissimilarity matrix
%
% DESCRIPTION
% A and B are either sets of cell structures of strings or strings.
% The edit distance is defined for strings of arbitrary lengths. It finds
% the minimum total cost of transforming one string into another, given
% the costs of insertion, deletion and substitution.
%
% Dynamic programming is used to compute the edit distance. Assume we are
% given two strings S and T of the lengths K and L, respectively. An
% (K+1)x(L+1) array M is filled with partial costs such that M(K,L)
% yields the edit distance D(S,T). The definition of the entries of M is
% recursive:
%   M(0,j) = j*C_del,  j=0,1,...,L
%   M(i,0) = i*C_ins,  i=0,1,...,K
% For the pairs (i,j), use
%   M(i,j) = min {M(i-1,j)+C_del, M(i,j-1)+C_ins, M(i-1,j-1)+C_subs(S(i)~=T(j)))
%
% DEFAULT
%   B    = A
%   COST = [1 1 1]
%   SC   = 1
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function d = editdistm(a,b,cost,scale)
if nargin < 4,
  scale = 0;
end
if nargin < 3,
  cost = [1 1 1];
end
if nargin < 2,
  b = a;
end

if length(cost) ~= 3,
  error('Wrong vector of edit costs.');
end
if any(cost) < 0,
  error('Edit costs should be nonnegative.');
end

if ~iscell(a),
  if isstr(a) | (~isstr(a) & min(size(a)) == 1)
    a = {a};
  else
    error('A is improper.');
  end
end
if ~iscell(b),
  if isstr(b) | (~isstr(b) & min(size(b)) == 1)
    b = {b};
  else
    error('B is improper.');
  end
end


m = length(a);
n = length(b);

d = zeros(m,n);
if scale == 1,
  for i=1:m
    al = length(a{i});
    for j=1:n
      bl = length(b{j});
      if bl < al, h = bl; bl = al; al = h; end
      scale = max([al*cost(3)+(bl-al)*cost(1), al*cost(2)+bl*cost(1),...
                   bl*cost(2)+al*cost(1), al*cost(3)+(bl-al)*cost(2)]);
      d(i,j) = editd(a{i},b{j},cost)/scale;
    end
  end
elseif scale == 0,
  for i=1:m
    for j=1:n
      d(i,j) = editd(a{i},b{j},cost);
    end
  end
else
  ;
end
return;




function d = editd (s,t,cost)
m = length(s);
n = length(t);

if isstr(s),
  s = ['a',s];
  t = ['a',t];
else
  if size(s,1) > size(s,2)
    s = [s(1); s];
    t = [s(1); t];
  else
    s = [s(1) s];
    t = [s(1) t];
  end
end

M = zeros(m+1,n+1);
M(1,1) = 0;
M(1,2:n+1) = (1:n) *cost(1);    % Insertions
M(2:m+1,1) = (1:m)'*cost(2);    % Deletions

for i=2:m+1
  for j=2:n+1
    M(i,j) = min ([M(i-1,j)+cost(2), M(i,j-1)+cost(1), M(i-1,j-1)+cost(3)*(s(i) ~= t(j))]);
  end
end
d = M(m+1,n+1);
