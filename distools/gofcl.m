%GOFCL  Goodness of clusters/classes separability vs compactness
%
%   J = GOFCL(D)
%
% INPUT
%   D	NxN Dissimilarity dataset
%
% OUTPUT
%   J	Criterion value
%
% DESCRIPTION
% Computes a goodness of clusters/classes in an NxN dissimilarity dataset D.
% D is labeled. The criterion provides a trade-off between the cluster/class 
% compactness and cluster/class separability. Consider K classes, with the total
% numbr of objects N and N_i elements in the i-th class. Let A_ij be the average 
% dissimilarity between the i-th and j-th clusters. Then the criterion is computed
% as:
%     J = sum_i (n_i sum_{j neq i} N_i/(N-N_i) A_ij)  / (2 sum_i n_i A_ii )
%
% The larger value, the better the separability between the classes.
%
% EXAMPLE
% Compare:
% 1) rand('seed',37); randn('seed',37); a=gendats(40,2,1); d=sqrt(distm(a)); 
%    scatterd(a); j1 = gofcl(d);
% 2) rand('seed',37); randn('seed',37); a=gendats(40,2,7); d=sqrt(distm(a)); 
%    scatterd(a); j2 = gofcl(d);


% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester

function J = gofcl(d);

issym(d);
lab = getnlab(d);
cc  = max(lab);

for i=1:cc
  Z = find(lab == i);
  dw(i) = mean(mean(d(Z,Z)));
	for j=i+1:cc
    ZZ = find(lab == j);
    db(i,j) =  mean(mean(d(Z,ZZ)));
    db(j,i) =  mean(mean(d(Z,ZZ)));
	end
end


J = 0;
for i=1:cc
  for j=i+1:cc
    J = J + db(i,j)/((dw(i)+dw(j)));
  end
end
J = J /((cc-1)*cc);
