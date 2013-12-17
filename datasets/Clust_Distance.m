function D = Clust_Distance(Xdat,type)
% Clust_Distance - Computes Distances between Samples
%
% <<< USE >>>
% D = Clust_Distance(Xdat,type)
%
% <<< INPUT >>>
% Xdat - (N x F) - Matrix of samples (One sample for each row)
% type - String  - Type of Distance Measure
%        {'Euclidean','PositiveCorrelation','MixedCorrelation','Angle'} or
%        {'e','p','m','a'}
%
% <<< OUTPUT >>>
% D - (N x N) - Matrix of distances. D(i,j) = distance between sample i and j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by E.P. van Someren  %%
%%  TU Delft           Jan 2000  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M] = size(Xdat);

switch type
 	case {'Euclidean','e'},							% EUCLIDEAN DISTANCE
       D=zeros(N,N);
       for i=1:N,
          D(i+1:end,i)=sqrt(sum((Xdat(i+1:end,:) - ones(N-i,1)*Xdat(i,:)).^2,2));
       end
       D=D+D';
 	case {'PositiveCorrelation','p'},			% PEARSON (POSITIVE) CORRELATION
		D = 1-corrcoef(Xdat');
	case {'MixedCorrelation','m'},				% MIXED CORRELATION
		D = 1 - abs(corrcoef(Xdat'));
	case {'Angle','a'},								% ANGLE BETWEEN VECTORS
		D = zeros(N);
		for i=1:N,
			for j = (i+1):N,
				a = Xdat(i,:)';
				b = Xdat(j,:)';
				D(i,j) = real((acos((a'*b)/(sqrt((a'*a))*sqrt((b'*b))))/(2*pi))*360);
				D(j,i) = D(i,j);
			end;
		end;
 	otherwise,
  		error(['EUSError: Unknown distance measure ' type '!']);
end;

% if the distance is zero, Matlab sometimes substitues -eps, which causes
% problems with log when we try to plot on the log scale.
neg_pos = find(D<eps);
D(neg_pos) = eps;
