function DB = hdb(a,linkage,dmeasure,maxNclusters)        

% -------------------------------------------------------
% DB = hdb(a,linkage,dmeasure,maxNclusters)        
%
% Apply Davies-Bouldin to a hierarchical clustering  
% 
% INPUT:
%       linkage: 'c' (complete);                      
%                's' (single);                        
%                'a'(average)                         
%      dmeasure: distance measure (See Clust\_Distance)  
%  maxNclusters: maximal number of clusters to evaluate 
%
% OUTPUT:
%   DB: Vector containing the Davies-Bouldin indeces
%       for the different numbers of clusters evaluated
%       (Note that the evaluation starts at two clusters
%       and runs up to maxNclusters)
% -------------------------------------------------------

x = +a;        
D = Clust_Distance(x,dmeasure);    
Xd = Clust_Dendrogram(D,linkage);    
for i = 2:maxNclusters,                
    l = Xd{3}(end-i+1,:);
    DB(i-1) = dbi(x,l);
end                                  
figure(100),clf                       
plot([2:1:maxNclusters],DB,'ro-')            
 