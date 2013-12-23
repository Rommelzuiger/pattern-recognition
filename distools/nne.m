%NNE Leave-one-out Nearest Neighbor Error on a Dissimilarity Matrix
%
%   [E,LAB] = NNE(D)
%
% INPUT
%   D   NxN symmetric dissimilarity dataset
%
% OUTPUT
%   E   Leave-one-out error
%   LAB Nearest neighbor labels
%
% DESCRIPTION
% Estimates the leave-one-out error of the 1-nearest neighbor rule
% on the givven symmetric dissimilairy data.
%

% Copyright: Robert P.W. Duin, r.p.w.duin@prtools.org and
% Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [e,NNlab] = nne(D)

[m,n] = size(D);
if m ~= n,
  error('Distance matrix should be square.');
end

lab = getlab(D);
[nlab,lablist] = renumlab(lab);
D(1:m+1:end)   = inf;
[d,M] = min(D');
e     = mean(nlab(M) ~= nlab);
NNlab = lablist(nlab(M),:);
return;
