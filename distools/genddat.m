%GENDDAT Generate random training and test sets for dissimilarity data
%
% 	[DTR,DTE,I,J] = GENDDAT(D,M,K)
%
% INPUT
% 	D 	NxN dissimilarity dataset
% 	M 	Cx1 vector of class sizes or frequencies, or a single number / frequency in(0,1)
% 	K 	Cx1 vector of class sizes or frequencies, or a single number / frequency in(0,1)
% 			(optional, default: K = M)
%
% OUTPUT
%   DTR Training dissimilarity dataset
%   DTE Test dissimilarity dataset
% 	I 	Indices of the training objects
% 	J 	Indices of the test objects
%
% DESCRIPTION
% Generates random training and test sets from a square dissimilarity dataset D.
% Feature labels and object labels of D should be equal. Note that M and K can be 
% either scalars or vectors with as many components as classes (=C), defining 
% specific sizes/fractions per class.
%
% Per default, all training are used as a represenation set (M=K). If M is Cx1 vector, 
% then GENDDAT selects at random M(i) vectors out of the i-th class in D and stores them 
% in the dataset DTR of the size [sum(M)]x[sum(M)]. The remaining  objects are stored in 
% [N-sum(M)]x[sum(M)] dissimilarity data DTE. Classes are ordered by using RENUMLAB(GETLAB(D)). 
% If M is a scalar, then M objects in total (given by number or frequency) are selected 
% at random according to the class priors. I and J are the indices of the training 
% and testing objects, respectively.
%
% If K is provided, then the first K(i) training objects per i-the class are used 
% for the representation set. Note that K(i) can be a frequency.
% 
% DEFAULT
% 	K = M
% 
% EXAMPLE
% Let D be 100 x 100 dataset with two classes [40 60] and class priors [0.4 0.6].
% 1)	[DTR,DTE] = GENDDAT(D,0.6) 
% DTR is 60x60 and DTE is 40x60. 60% of objects of the first class (24 in total) 
% and 40% of the objects of the second class (36 in total) are selected for DTR.
%
% 2) [DTR,DTE] = GENDDAT(D,0.6,0.1) 
% DTR is 60x6 and DTE is 40x6. 60% of objects of the first class (24 in total) 
% and 40% of the objects of the second class (36 in total) are selected for training.
% From that, 10% of the first training objects per class are selected for the 
% represenatation set. 10% from 24 rounds to 2 objects for the first class, while
% 10% of 36 rounds to 4 objects for the second class. 
%
% SEE ALSO
% DATASETS, RENUMLAB
%

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org, and
% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology, and
% School of Computer Science, University of Manchester


function [DTR,DTE,Itr,Ite] = genddat(D,m,k);

[n,nk,c] = getsize(D);
nlab     = getnlab(D);
discheck(D,[],1);   % allow for similarities     

if nargin < 3, 
	k = []; 
else
	if length(k) == 1
		k = k*ones(1,c);
	elseif length(k) == c
		;
	else
		error('Vector length of the number of objects should equal the number of classes.')
	end
	
	if ~(all(k == round(k))) & ~(all(k > 0 & k < 1))
		error('K should be given either by integers or frequencies in (0,1).')	
	end
end  
if ~(all(m == round(m))) & ~(all(m > 0 & m < 1))
	error('M should be given either by integers or frequencies in (0,1).')	
end

[ja,jb] = gendat(dataset([1:n]',nlab),m);
ja = +ja; 
jb = +jb;

J  = [];
for j=1:c
	K = find(nlab(ja)==j);
	if ~isempty(k)
		if k(j) < 1, 
			k(j) = round(k(j)*length(K)); 
		end
		if k(j) > length(K)
			error('Requested size of the representation set is not possible.')
		end
		K = K(1:k(j));
	end
	J = [J; K(:)];
end
M = setdiff([1:length(ja)]',J);

Itr = [ja(J); ja(M)];
Ite = jb;
DTR = D(Itr,ja(J));
DTE = D(Ite,ja(J));

