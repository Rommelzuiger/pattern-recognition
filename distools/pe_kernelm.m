%PE_KERNELM Pseudo-Euclidean kernel mapping
% 
%    K = PE_KERNELM(A,B)
%    W = B*PE_KERNELM
%    W = PE_KERNELM([],B)
%    K = A*W
%
% INPUT
%   A     Pseudo-Euclidean dataset of size NxK
%   B     Pseudo-Euclidean dataset of size MxK
%
% OUTPUT
%   W     PE mapping
%   K     Kernel matrix, size [N M]
%
% DESCRIPTION
% Computation of a kernel matrix in a pseudo-Euclidean space. The signature
% of this space should be stored in the datasets A and B, see SETSIG.
% K = A*J*B', where J is a diagonal matrix with 1's, followed by -1's.
% J = diag ([ONES(SIG(1),1);  -ONES(SIG(2),1)]);
% The two-element vector SIG stores the signature of the space. This is the
% number of 'positive' dimensions, followed by the number of 'negative'
% dimensions. It is computed by a pseudo-Eucledean embedding, e.g. PSEM,
% and stored in the related mapping and datasets that are projected in this
% space.
%
% The resulting kernel matrix K is indefinite in case A == B. This routine
% may be used in support vector routines and other kernelized procedures.
% Note that most of such routines are not optimal for indefinite kernels.
%
% EXAMPLE
% a = gendatb;                  % generate dataset
% d = a*proxm(a,'m',1);         % compute L1 distance matrix
% w = psem(d);                  % embed in PE space
% b = d*w;                      % project data in this space
% [trainset testset] = gendat(b,0.5);     % split in trainset and testset
% ktrain = pe_kernelm(trainset,trainset); % compute train kernel
% w = svc(ktrain,0);            % compute SV classifier
% ktest = pe_kernelm(testset,trainset);   % compute test kernel
% ktest*w*testc                 % inspect error of testset   
%
% SEE ALSO
% DATASETS, MAPPINGS, PE_EM, PE_DISTM

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = pe_kernelm(a,b);

  if nargin < 2, b = []; end
  mapname = 'PE kernel mapping';
	if (nargin < 1) | (isempty(a))  
		% Definition: pe kernel mapping.
		w = mapping(mfilename,{b});
		w = setname(w,mapname);
	
  elseif isdataset(a)
    if isempty(b)   % store a as rep set, 'training'
      w = mapping(mfilename,'trained',a,getlab(a),size(a,2),size(a,1));
		  w = setname(w,mapname);
    elseif isdataset(b); % compute kernel between a and b
      w = pe_mtimes(a,b');
    elseif ismapping(b)
      if isuntrained(b) % nothing stored yet, do it now, a is rep set
        w = mapping(mfilename,'trained',a,getlab(a),size(a,2),size(a,1));
		    w = setname(w,mapname);
      else % we have already a rep set: compute kernel matrix
        b = getdata(b,1);
        w = pe_mtimes(a,b');
      end
    else
      error('Second parameter should be dataset')
    end
    
  else  % may be double ???
    a = dataset(a);
    w = feval(mfilename,a,b);
    
  end
  
return
      
      