function [P,lab,WS] = kmclust(a,k,plotflag,stepflag)

%--------------------------------------------------
% [P,lab,WS] = kmclust(a,k,plotflag,stepflag)
%
% KM Performs K-means clustering
% 
% INPUT: 
%          a: (n X p) data set
%          k: number of desired clusters
%   plotflag: if set plots all the samples,
%             the final positions of the 
%             prototypes as well as the 
%             trajectories of the prototypes.
%   stepflag: if set, the function single steps
%             through the k-means algorithm, 
%             displaying the positions of the 
%             prototypes and the samples, as
%             well as the assignment of the 
%             samples to the prototypes.
% OUTPUT:
%          P: (k X m) matrix, the prototypes
%       labs: (n X 1) vector, the cluster labels
%         WS: the within-scatter associated with
%             the obtained clustering.
%--------------------------------------------------

X = +a;
thresh = 1e-6;
[n,m] = size(X);
cols = {'r','g','k','m','c','b','r','g','k','m','c','b','r','g','k','m','c','b','r','g','k','m','c','b'};

% STEP 1: GENERATE THE PROTOTYPES: RANDOMLY PICK k SAMPLES FROM DATA SET
r = randperm(n);
P = X(r(1:k),:);

%Mx = max(X,[],1);
%Mn = min(X,[],1);
%P = rand(k,m).*(ones(k,1)*(Mx - Mn)) + ones(k,1)*Mn;

Ppos(:,:,1) = P;      % positions of the prototypes
ii = 1;

WS = [];
while 1,
   
   ii = ii + 1;
   %dprot = distmat(X,P);
   dprot = distm(a,P);

   [min_dprot,lab]=min(+dprot,[],2);
   
   if plotflag,
      figure(5),clf
      %plot(X(:,1),X(:,2),'+')
      hold on
      for i = 1:k,
         c = X(find(lab == i),:);
         plot(c(:,1),c(:,2),['.',cols{i}],'markersize',25)
      end
      plot(P(:,1),P(:,2),'+','markersize',20,'linewidth',5)
   end
   
   % STEP 3: COMPUTE NEW PROTOTYPES
   Pnew = [];
   for i = 1:k,
      samples_in_clust = X(find(lab == i),:);
      Pnew(i,:) = mean(samples_in_clust,1);
   end
   
   % STEP 4: CHECK CONVERGENCE
   Pchange = mean((P - Pnew).^2,2);
   maxPchange = max(Pchange);
   if maxPchange < thresh,
      break
   end
   
   % STEP 5: UPDATE PROTOTYPES
   P = Pnew;
   
   Ppos(:,:,ii) = P;
   if plotflag & stepflag,
       pause
   end
end

% WITHIN SCATTER
WS = sum(min_dprot);

if plotflag,
   for i = 1:k,
      for j = 1:size(Ppos,3),
         xtraj(j) = Ppos(i,1,j);
         ytraj(j) = Ppos(i,2,j);
      end
      plot(xtraj,ytraj,'k.-','linewidth',2)
   end
end