%CHECKEUCL Check whether a square dissimilarity matrix has a Euclidean behavior
%
%   [NEF,NER,W] = CHECKEUCL(D)
%   [NEF,NER]   = CHECKEUCL(D,K)
%   [NEF,NER,K] = CHECKEUCL(D,'all')
%
% INPUT
%   D     NxN dissimilarity matrix or dataset
%   K     Vector with desired subset sizes
%
% OUTPUT
%   NEF   Index of non-Euclidean behavior; negative eigen-fraction NEF in [0,1)
%   NER   Index of non-Euclidean behavior; negative eigen-ratio NER >= 0
%   W     Pseudo-Euclidean embedding
%
% DESCRIPTION
% Computes how well the square dissimilarity matrix D can be embedded
% in a Euclidean space.
% NER and NEF are computed by performing the Pseudo-Euclidean embedding.
% NEF is a ratio expressing the sum of magnitudes of all negative
% eigenvalues divided by the sum of magnitudes of all eigenvalues. NER is
% a ratio of the magnitude of the lowest negative eigenvalue to the largest
% positive eigenvalue.
%
% D is Euclidean if D.^2 is isometrically embeddable into a Euclidean space.
% Ideally, both NEF and NER are zero. Note that due to numerical inaccuracies
% of the emebdding procedure, both NEF and NER might be very small, e.g.
% in the order of ~1e-10, for perfect Euclidean distance data.
%
% In case a set of subset sizes K is given as many random subsets of K(i) 
% objects are selected that are needed to estimate the expected NEF for
% matrices of K(i)*K(i) dissimilarities with a standard deviation of 5%.
%
% SEE ALSO
% PE_EM

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester
%

function [nef, ner, W] = checkeucl(D,N)

	if nargin < 2, N = []; end
	alf = 0.05;
	[m,k] = size(D);
	if m ~= k
		error('Dissimilarity matrix D should be square.')
	end

	D = dataset(D,1);  % we are not interested in labels here.
	D = setfeatlab(D,ones(m,1));
	
	if isempty(N)
		W    = pe_em(D);
    L    = getdata(W,'eval');
		nef  = sum(abs(L(find(L < 0))))/sum(abs(L));
		ner  = max(abs(L(find(L < 0)))) / max(L(find(L > 0)));
    if isempty(ner), ner = 0; end
  else
    if isstr(N) & strcmp(N,'all')
      N = [3,4,5,7,10,15,20,30,50,70,100,150,200,300,500,700,1000,1500,2000,3000];
      N = [N(N<m) m];
    end
    npoints = length(N);
		nef = zeros(1,npoints);
		ner = zeros(1,npoints);
    if npoints > 1
      s = sprintf('checkeucl: NEF on %i points: ',npoints);
      prwaitbar(npoints,s);
    end
		for j = 1:npoints
      if npoints > 1, prwaitbar(npoints,j,[s int2str(j)]); end
			n = N(j);
      if n==m
        [nef(j) ner(j)] = feval(mfilename,D);
      else
        nfe = 0;
        nre = 0;
        nfv  = 0;
        for i=1:5
          [nf nr] = feval(mfilename,genddat(D,n));
          nfe = nfe+nf;    nfee = nfe/i;
          nre = nre+nr;    
          nfv = nfv+nf*nf; nfvv = nfv/i;
        end
        acc = sqrt(nfvv-nfee^2)/(sqrt(i)*nfee);
        while acc > alf & i < 1000
          i = i+1;
          [nf nr] = feval(mfilename,genddat(D,n));
          nfe = nfe+nf;    nfee = nfe/i;
          nre = nre+nr;    
          nfv = nfv+nf*nf; nfvv = nfv/i;
          acc = sqrt(nfvv-nfee^2)/(sqrt(i)*nfee);
          if nfee < 0.001, nfee = 0; acc = 0; end
          %disp([i,round(10000*acc)])
        end
        nef(j) = nfee;
        ner(j) = nre/i;
      end
		end
    if npoints > 1, prwaitbar(0); end
    W = N;
	end

return;
