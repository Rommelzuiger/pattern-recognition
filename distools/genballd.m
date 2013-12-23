%GENBALLD Generate ball distance matrix
%
%   D = GENBALLD(N,K,E)
%
%  N 1xC vector with number of objects per class
%  K Dimensionality
%  E 1xC with ball sizes, default: class number

function d = genballd(n,k,e)

c = length(n); % number of classes
if nargin < 2, k = 1; end
if nargin < 3, e = [1:c]; end
m = sum(n);

x = 100*rand(m,k);
r = [];
for j=1:c
    r = [r; e(j).*ones(n(j),1)];
end

d = sqrt(distm(x)) - repmat(r,1,m) - repmat(r',m,1); 
d = d - diag(diag(d));
attempts = 0;
while(any(d(:)<0))
for j=2:m
    while any(d(j,:) < 0)
				attempts = attempts + 1;
				if attempts > 50*m
					error('No solution found, shrink ball sizes or number of balls, or enlarge dimensionality')
				end
        x(j,:) = 100*rand(1,k);
        d(j,:) = sqrt(distm(x(j,:),x)) - repmat(r(j),1,m) - repmat(r',1,1);
        d(:,j) = d(j,:);
        d(j,j) = 0;
    end
end

end         

labs = genlab(n);
d = dataset(d,labs);
d = setfeatlab(d,labs);
d = setprior(d,0);