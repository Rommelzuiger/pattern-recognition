%DISSTAT Basic information (statistics) on a dissimilarity represenatation 
%
%   [MIN,MAX,ME,STD,SKEW] = DISSTAT(D)
%
% INPUT
%   D   NxN or NxM dissimilarity matrix or dataset
%
% OUTPUT
%   MIN  Vector of minimal dissimilarities per class and for the whole data   
%   MAX  Vector of maximal dissimilarities per class and for the whole data   
%   ME   Vector of average dissimilarities per class and for the whole data   
%   STD  Vector of standard deviations of dissimilarities per class and for the whole data  
%   SKEW Vector of skeweness coefficients per class and for the whole data
%
% DESCRIPTION
% Basic statistics on a dissimilarity represenatation: minimum, maximum,
% average, standard deviation and skewness are computed per class and for
% the whole data. They are returned in vectors of the length of C+1, where C
% is the number of classes. The last elements reflect the statistics for
% the whole data.

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com, Robert P.W. Duin, r.duin@ieee.org
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [MIN,MAX,ME,STD,SKEW] = disstat(D)

if isnumeric(D)
    DD   = D;
    DD   = setdiff(DD(:),0);
    MIN  = min(DD);
    MAX  = max(DD);          
    ME   = mean(DD);
    STD  = std(DD);
    SKEW = mean(((DD-ME)/STD).^3);
else    
    ll  = getlab(D);
    fll = getfeatlab(D);
    [lab,flab,lablist] = renumlab(ll,fll);

    D = +D;
    I = find(D == 0);

    for i=1:max(lab)+1
        I = find(lab == i);
        J = find(flab == i);
        if i < max(lab)+1
            DD = D(I,J);
        else
            DD = D;
        end
        DD     = setdiff(DD(:),0);
        MIN(i,1) = min(DD);
        MAX(i,1) = max(DD);          
        ME(i,1)  = mean(DD);
        STD(i,1) = std(DD);
        SKEW(i,1)= mean(((DD-ME(i))/STD(i)).^3);
    end
end
return;