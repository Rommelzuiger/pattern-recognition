%STRKERM  String Kernel Matrix by Lodhi et al
%
% 	[K,KK] = STRKERM (A,B,W,LAM,KNORM)
%
% INPUT
%   A     Cell structure of N strings or a single string
%   B     Cell structure of M strings or a single string (optional; default: B = A)
%   W     Vector of K weights (optional; default: [1 1])
%   LAM   Decay factor Weight parameter in (0,1] used to scale the intermediate kernels
%         (optional; default: 1)
%   KNORM Parameter (0/1) indicating whether the kernel should be normalized
%         to remove bias wrt to text length (optional; default: 1)
%
% OUTPUT
%   K 		NxM kernel matrix
%   KK		NxMxK matrix of intermediate kernels
%
% DESCRIPTION
% Derives a string kernel matrix K following the idea of Lodhi et al [2002]. A and B are 
% either individual strings or cell structures of strings. The basic idea is to compute 
% the kernel as an inner product of the feature vectors for two strings and give a sum 
% over all common subsequences weighted according to their frequency of occurrence and 
% lengths. Dynamic programming is used to compute the kernel values between strings. 
% If dataset is needed, it should be created afterwards.
%
% DEFAULT
%   B     = A
%   W     = [1 1]
%   LAM   = 1
%   KNORM = 1
%
% REFERENCE
% H.Lodhi, C.Saunders, J. Shawe-Taylor, N.Cristianini, C.Watkins, "Text
% Classification using String Kernels", J. of Machine Learning Research 2,
% 419-444, 2002.
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% EWI Faculty, Delft University of Technology and
% School of Computer Science, University of Manchester


function [K,KK] = strkerm (a,b,w,lam,Knorm)
if nargin < 5,
  Knorm = 1;
end
if nargin < 4 | isempty(lam),
  lam = 1;
end
if nargin < 3 | isempty(w),
  w = [1 1];
end
if nargin < 2 | isempty(b),
  b = a;
end

if any(w) < 0,
  error('Weights should be nonnegative.');
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


K = zeros(m,n);
for i=1:m
  for j=1:n
    K(i,j) = strkernel(a{i},b{j},w,lam);
  end
end

if Knorm,
  for i=1:m
    Kaa(i,1) = strkernel(a{i},a{i},w,lam);
  end
  for j=1:n
    Kbb(1,j) = strkernel(b{j},b{j},w,lam);
  end
  K = K ./ sqrt(Kaa*Kbb);
end
return;






function K = strkernel(s,t,w,lam)

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

ns = length(s);
nt = length(t);
n  = length(w);

K = 0;
KK(:,:,1) = ones(ns+1,nt+1);
for i=1:n
  KK(:,:,i+11) = zeros(ns+1,nt+1);
  for j=1:ns
    ss = 0;
    for k=1:nt
      if t(k) == s(j),
        ss = ss + KK(j,k,i);
      end
      KK(j+1,k+1,i+1) = KK(j,k+1,i+1) +ss;
    end
  end
  K = K + w(i)*KK(ns+1,nt+1,i+1)*lam^(2*i);
end
