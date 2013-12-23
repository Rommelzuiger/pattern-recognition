%KCENTERM Kernel Weighted Centering Mapping
%
%   W = kcenterm(K,V)
%
% INPUT
%   K   NxN kernel or similarity matrix (dataset)
%   V   'C' (default) or Nx1 weight vector of nonnegative values; sum(V)=1
%
% OUTPUT
%   W   Kernel weighted centering mapping
%
% DEFAULT
%   V   = 'C';
%
% DESCRIPTION
% Defines a mapping that centers the given NxN kernel matrix K:=K(X,X)
% such that the weighted mean coincides with the origin in the vector space
% induced by K. The nonnegative weights of the mean are defined by an
% Nx1 vector V; sum(V) = 1. V = 'C' stands for V=ones(N,1)/N, hence
% the true centering.
%
% The mapping works with both positive (semi-)definite and indefinite
% kernels. V can be applied to a new MxN kernel Knew:=K(Xnew,X).
%
% SEE ALSO
% KPCA, MAPPINGS, DATASETS
%

% Copyright: Ela Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester
%

function [V,KK] = kcenterm(K,v)

if nargin < 2 | isempty(v)
  v = 'c';
end

if nargin < 1 | isempty(K)
   V = mapping(mfilename);
   V = setname(V,'kcenterm');
   return
end


if (isdataset(K) | isa(K,'double'))
  if ismapping(v)

    [m,n] = size(K);
    pars  = getdata(v);

    v   = pars{1};
    Kwm = pars{2};

    H  = eye(n) - v*ones(1,n);
    V  = (K - repmat(Kwm',m,1))* H;
    return;
  end
end



lab      = getlab(K);
lablist  = getlablist(K);
[n,m,c]  = getsize(K);

tol = 1e-12;
if ~issym(K,tol),
  error('Kernel matrix K should be symmetric.')
end

if v == 'c',
  v = ones(n,1)/n;
end

if length(v)==1
  l = intersect(v,1:n)
  v = zeros(n,1);
  v(l) = 1;
end

if length(v) ~= n,
  error('V has a wrong size.');
end

if any(v) < 0 | any(v) > 1,
  error('V should have elements in [0,1].');
end

if abs(sum(v) - 1) > tol
  error('sum(V) ~= 1.');
end


% Center K such that the weighted mean coincides with the origin
Kwm  = K*v;
if nargout > 1
  H  = eye(n) - v*ones(1,n);
  KK = H * K * H;
end

V = mapping(mfilename,'trained',{v,Kwm},[],m,m);
V = setname(V,'kcenterm');
return
