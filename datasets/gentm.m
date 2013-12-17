function gentm(a,nnodes,nbasis,reg,niter)

% GENTM(A,NNODES,NBASIS,REG,NITER)
%
% Wrapper function for the generative 
% topographic map (GTM) (See GTM for more
% detail.
%
%      A: dataset
% NNODES: vector with number of nodes per dimension
%         e.g. [5 7]. If a scalar is given as input,
%         the same number of nodes is employed in the 
%         other dimension(s)
% NBASIS: vector with number of basis functions per
%         dimension, e.g. [3 5]. If a scalar is given 
%         as input, the same number of nodes is employed 
%         in the other dimension(s) 
%    REG: regularization parameter
%  NITER: number of iterations

if min(min(nbasis))<2,
    error('At least two basis function per dimension required')
end
if min(min(nbasis))<2,
    error('At least two nodes per dimension required')
end
n = size(getdata(a),1);
b = a-ones(n,1)*mean(a);
for i = 1:niter
  w = gtm(b(:,1:2),nnodes,nbasis,'mean',reg,1e-5,i); 
  figure(1), clf; scatterd(b); plotgtm(w); drawnow; 
  pause(0)
end;

