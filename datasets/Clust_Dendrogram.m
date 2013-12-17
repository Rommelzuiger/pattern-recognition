function Xd = Clust_Dendogram_new(D,type)

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% <<<<<<<<<<<<< Version employed to update Average Linkage >>>>>>>>>>>>>>>
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% Clust_Dendogram - Computes Dendogram given Distance Matrix
%
% <<< USE >>>
% Xd = Clust_Dendogram(D,type)
%
% <<< INPUT >>>
% D    - (N x N) - Distance Matrix (D(i,j) = distance between sample i and j)
% type - String  - Type of Linkage
%        {'Single','Average','Complete'} or
%        {'s','a','c'}
%
% <<< OUTPUT >>>
% Xd{1} = Cl    - correct object order for plotting
% Xd{2} = Conn  - connection history for plotting
% Xd{3} = CHist - complete record of clustering process

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by L.F.A.Wessels           %%
%%  TU Delft      November 1999        %%
%%  Modifications by P.J. van der Veen %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% size
[Nobjects, Nclusters] = size(D);

% if not square, error
if (Nobjects ~= Nclusters),
 error('distance matrix not square')
end

% put NaNs on the diagonal: distance from object to itselft irrelevant
D = D + diag(NaN * ones(Nobjects, 1), 0);

% save cluster ORDER information in Cl
for i = 1:Nobjects, Cl{i} = uint16(i); end
if Nobjects > 65535
 error('Nobjects > 65535');
end

% Cluster history: all partitions are remebered
CHist = zeros(Nobjects,Nobjects);
CHist(1,:) = 1:Nobjects;

% count of the number of objects in a cluster, for average linkage
Z = ones(1, Nobjects);

% find objects closest to each other
% find min over columns and save row indeces where this occurs
[min_over_columns, min_row_indeces] = min(D, [], 1);

% store the available elements
ellist = 1:Nobjects;

%---COMPUTE DENDROGRAM------------------------------------------------
% loop until Nclusters == 1
for K = 1:Nobjects - 1,

 if (mod(K,100)==0)
  disp([' Done ',num2str(K),' / ',num2str(Nobjects)])
 end

 % find min of the column minima; the index is the column pos. of
 % the matrix minimum
 [the_min, min_col_index] = min(min_over_columns);

 % the column minimum index is used to determine the row position
 % of the matrix minimum
 min_row_index = min_row_indeces(min_col_index);
 
 % rank
 I = min(min_row_index, min_col_index);
 J = max(min_row_index, min_col_index);

 % rank (recalculated)
 I2 = ellist(I);
 J2 = ellist(J);

 % Connection history
 % save the connection information in Conn
 % at each iteration, the 1st objects of the
 % clusters being merged are listed. Eg if
 % at iteration i, clusters 2 and 3 are merged and
 % the first elements of clusters 2 and 3 are 5 and
 % 6 respectively, with the distance between them D, then
 % Conn(i,:) = [5, 6, D, 2, 3]
 Conn(K,:) = [double([Cl{I}(1), Cl{J}(1)]), the_min, double([I2, J2])];
 
 % here we know that clusters/objects I and J need to be merged
 Cl{I} = [Cl{I}, Cl{J}];

 % Cluster history: all partitions are remebered
 CHist(K+1,:) = CHist(K,:);
 CHist(K+1,find(CHist(K+1,:)==J))=I;

 % update number of clusters
 Nclusters = Nclusters - 1;

 % recompute distances between new cluster and the rest
 switch type, 
  case {'Complete','c'},					% Complete 
   D(:,I) = max(D(:,I), D(:,J));
   D(I,:) = max(D(I,:), D(J,:));
   D(I,I) = NaN;
  case {'Single','s'},						% Single 
   D(:,I) = min(D(:,I), D(:,J));
   D(I,:) = min(D(I,:), D(J,:));
   D(I,I) = NaN;
  case {'Average','a'},						% Average
   % average of distances between clusters
   D(I,:) = (Z(I) * D(I,:) + Z(J) * D(J,:)) / (Z(I) + Z(J));
   % update count of number of objects in a cluster
   Z(I) = Z(I) + Z(J); Z(J) = 0;
  otherwise 
   error('Unknown linkage type')
 end
 
 %disp(['merged: ',num2str(i),' ',num2str(J)])
 %printmat(D)
 
 % Update the distance matrix resulting from combining the closest 
 % objects: remove the row and column J
 D(:, J) = inf; D(J, :) = inf;

 ellist(J+1:end) = ellist(J+1:end) - 1;  
 ellist(J) = 0;

 % recompute the mean_over_columns
 el = find(min_row_indeces == J | min_row_indeces == I);
 [min_over_columns(el), min_row_indeces(el)] = min(D(:,el), [], 1);
 min_over_columns(J) = inf; 
 
end % for K = 1:Nobjects - 1

% return the correct object order
Cl = Cl{1}; 
% pack variables in one cell struct
Xd{1} = Cl;
Xd{2} = Conn;
Xd{3} = CHist;


