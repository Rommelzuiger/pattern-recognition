% LAPLACE  Laplacian distributed random numbers.
%    LAPLACE(N) is an N-by-N matrix with random entries, chosen from a
%    Laplacian distribution with mean zero and variance one.
%    LAPLACE(N,M) is an N-by-M matrix with random entries.
%    LAPLACE(N,M,MU,S) is an N-by-M matrix with random entries, chosen 
%    from a Laplacian distribution with mean MU and covariance matrix S.
%    MU should be an 1-by-M vector, S an M-by-M matrix.
%    LAPLACE with no arguments is a scalar whose value changes each time
%    it is referenced. 

function out = laplace (n, m, mu, S)

	if (nargin < 1), n  = 1;          end;
	if (nargin < 2), m  = n;          end;
	if (nargin < 3), mu = zeros(1,m); end;
  if (nargin < 4), S  = eye(m);     end;

	out = myexprnd(1,n,m); 

	% Convert exponential to Laplacian distributed data

	for i = 1:n
		for j = 1:m
			if (rand(1,1) > 0.5)
				out(i,j) = -out(i,j);
			end;
		end;
	end;

	% Remove covariance
	out = out * inv(sqrtm(cov(out)));

	% Add in desired covariance and mean
	out = out * sqrtm(S) + ones(n,1)*mu;

return
