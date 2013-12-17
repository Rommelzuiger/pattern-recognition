function DB = dbi(x,l)

%-----------------------------------
% DB = dbi(x,l)
%
% Computes the Davies-Boulden index
% for the clustering defined by the 
% labels in l for the objects in the 
% data matrix x.
%
% INPUT: 
%  x: (n x p) data matrix
%  l: (n x 1) cluster labels
%
% OUTPUT:
%  DB: Davies-Bouldin index for the
%      defined clustering 
%-----------------------------------


ul = unique(l);   % set of cluster labels
c = length(ul);   % number of clusters
n = [];
for i = 1:length(ul),
    n(i) = sum(l==ul(i));
end

%---compute cluster means 
m = [];
for i=1:c,
    m(i,:) = mean(x(l==ul(i),:),1);
end
m;

%---compute cluster within scatter (scalars)
s = [];
for i=1:c,
    s(i) = sqrt(mean(sum((x(l==ul(i),:) - ones(n(i),1)*m(i,:)).^2,2)));
end
s;

%---compute the R-matrix
R = zeros(c,c);
for i = 1:c,
    for j = i+1:c,
        R(i,j) = (s(i) + s(j))/(sqrt(sum((m(i,:) - m(j,:)).^2)));
    end
end
R = R + R';   % same upper and lower diagonal
R;

%---most similar cluster
Rm = max(R,[],2);

%---compute DB
DB = mean(Rm);

return
