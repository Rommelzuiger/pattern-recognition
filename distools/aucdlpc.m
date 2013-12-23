%AUCDLPC AUC-LP classifier on dissimilarity data
%
%   [W,R] = AUCLPC(D,C)
%
% INPUT
%   D   Dissimilarity data
%   C   Trade-off parameter in the LP-optimization, similar to C in the SVM
%       (optional; default: 10)
%
% OUTPUT
%   W   Mapping: AUC-LP Classifier
%   R   Object identifiers of support objects
%
% DEFAULTS
%   C  = 10
%
% DESCRIPTION
% Determines a linear programming classifier for the dissimilarity representation D.
% This LP-function is determined for a two-class problem based on the maximization
% of the AUC measure (area under the curve) of the corresponding ROC curve.
% It calls AUCLPM(D,C,'SUBK',1).
%
%
% SEE ALSO
% MAPPINGS, DATASETS, AUCLPC, AUCLPM
%
% REFERENCE
% D.M.J. Tax and C. Veenman, Tuning the hyperparameter of an AUC-optimized classifier,
% in: Seventeenth Belgium-Netherlands conference on artificial intelligence, 2005, 224-231.

% David Tax, Robert P.W. Duin, Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty of Electrical Engineering, Mathematics and Computer Science,
% Delft University of Technology, The Netherlands.


function [w,r] = aucdlpc(d,C,usematlab,prec)

if nargin < 4, prec = 1e-7; end;
if nargin < 3, usematlab = 0; end
if nargin < 2 | isempty(C),
  C = 10;
end

if nargin < 1 | isempty(d)
    w = mapping(mfilename,C);
    w = setname(w,'AUCDLPC');
    return
end


[m,k,c] = getsize(d);

if c > 2,
  % multi-class problem

  lab = getlab(d);
  fe  = getfeatlab(d);
  [nlab,nfe,lablist] = renumlab(lab,fe);

  w = [];
  N = [];
  for i=1:c
    mlab = 2 - (nlab == i);
    mfe  = 2 - (nfe == i);
    dd   = setlabels(d,mlab);
    dd   = setfeatlab(dd,mfe);
    if ~isempty(d.prior)
      dd = setprior(dd,[d.prior(i),1-d.prior(i)]');
    end
    [v,J]= aucdlpc(dd,C,usematlab,prec);
    w = [w,setlabels(v(:,1),lablist(i,:))];
    N = [N; J];
  end
  r = unique(N);
  w  = setname(w,'AUCDLPC');
  return

else
  % two-class problem

  v    = auclpm(d,C,'subk',1,1,usematlab);
  nlab = getnlab(d);
  p    = getprior(d);
  q    = +(d*v);
  [qq,J] = sort(q);
  n1 = sum(nlab == 1);
  n2 = sum(nlab == 2);
%  e = n1-n2;
%  J1 = cumsum(nlab(J) == 1);
%  J2 = n2-cumsum(nlab(J) == 2);
  e1 = p(1)*cumsum(nlab(J) == 1)/n1 + p(2)*(1 - cumsum(nlab(J) == 2)/n2);
  e2 = p(2)*cumsum(nlab(J) == 2)/n2 + p(1)*(1 - cumsum(nlab(J) == 1)/n1);

  [min1,k1] = min(e1);
  [min2,k2] = min(e2);
% if min1 < min2          % NO!! We should assume that AUCLPM gives a proper
%    if k1 == length(q)   % and well defined class order: 1 -> -, 2 > +
%      w0 = q(J(k1)) + realmin;
%    else
%      w0 = (q(J(k1))+q(J(k1+1)))/2;
%    end
%    w = v.data.u - v.data.v;
%  else
    if k2 == length(q)
      w0 = q(J(k2)) + realmin;
    else
      w0 = (q(J(k2))+q(J(k2+1)))/2;
    end
    w = v.data.v - v.data.u;
%  end

  ss = sum(abs(w));
  J  = find(abs(w) >= ss*prec);
  w  = featsel(k,J)*affine(w(J),w0,d(:,J),getlablist(d),k,c);
  w  = setout_conv(w,1);
  w  = setname(w,'AUCDLPC');
  r  = J;
end
return
