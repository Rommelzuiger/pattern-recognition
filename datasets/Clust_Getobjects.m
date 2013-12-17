function [o2c,c2o] = Clust_Getobjects(Nc,Xd)


% --------------------------------------------------------------------
% o2c = getobjects(Nc,CHist)
%
% Returns the clustering for a given number of clusters
%
% INPUT:
%  Nc: number of clusters
%  CHist: clustering history in the new format:
%         (Nobjects + 1) X (Nobjects)
%         column i represents object i, and 
%         each row the cluster it is associated with.
%
% OUTPUT:
%  o2c: colunm vector of length Nobjects containing
%       the cluster index for object i at position i
% --------------------------------------------------------------------

% CHist is a Nobject X Nobjects matrix, with each row containing,
% at column i, a cluster label for object i. Note that the set of
% cluster labels at a given row are not necessarily a consecutive
% set of numbers.

CHist = Xd{3};

% number of objects clustered
Nobjects = size(CHist,2);
% error check
if Nc < 1,
 error('Nc should be larger than 0')
end
if Nc > Nobjects,
 error(['Nc should be smaller than Nobjects (',(num2str(Nobjects)),')'])
end
% row in CHist associated wit Nc clusters
row = CHist(Nobjects - Nc + 1, :);
% set of unique cluster labels at this level
U = sort(unique(row));
% create o2c
o2c = zeros(Nobjects,1);
for i = 1:length(U),
 o2c(find(row == U(i))) = i;
end

% generate c2o (cluster to objects)
c2o = [];
uniq = unique(o2c);
for i=1:length(uniq),
 c2o{i} = find(o2c == uniq(i));
end


