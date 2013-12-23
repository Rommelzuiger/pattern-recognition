%PROXXM Proximity Mapping
%
%   W = PROXXM(A,TYPE,P,WEIGHTS)
%   W = A*PROXXM([],TYPE,P,WEIGHTS)
%
% INPUT
%   A       MxK Dataset
%   TYPE    Type of the proximity (optional; default: 'distance')
%   P       Parameter of the proximity (optional; default: 1)
%   WEIGHTS Weights (optional; default: all 1)
%
%   OUTPUT
%   W       Proximity mapping
%
% DESCRIPTION
% Computation of the KxM proximity mapping (or kernel) defined by the
% MxK dataset A. Unlabeled objects in A are neglected. If B is an NxK
% dataset, then D=B*W is the NxM proximity matrix between B and A. The
% proximities can be defined by the following types:
%
%   'POLYNOMIAL'   | 'P':  SIGN(A*B').*(A*B').^P
%   'HOMOGENEOUS'  | 'H':  SIGN(A*B').*(A*B').^P
%   'EXP'          | 'E':  EXP(-(||A-B||)/P)
%   'EXP-SUM'      | 'ES': SUM_Z SIGN(P(Z)) * EXP(-(||A-B||)/P(Z))
%   'RBF'          | 'R':  EXP(-(||A-B||.^2)/(P*P))
%   'RBF-SUM'      | 'RS': SUM_Z SIGN(P(Z)) * EXP(-(||A-B||.^2)/(P(Z)^2))
%   'SIGMOID'      | 'S':  SIGM(A*B'/P)
%   'DSIGMOID'     | 'DS': SIGM(||A-B||^(2P(1))/P(2))
%   'DISTANCE'     | 'D':  ||A-B||.^P
%   'MULTIQUADRIC' | 'MQ': SQRT(||A-B||.^2/P(1) + P(2))
%   'THIN-PLATE'   | 'TP': ||A-B||.^(2*P(1))/P(2) * LOG(1+||A-B||.^(2*P(1))/P(2))
%   'N-THIN-PLATE' | 'NTP': -||A-B||.^(2*P(1))/P(2) * LOG(1+||A-B||.^(2*P(1))/P(2))
%   'MINKOWSKI'    | 'M':  SUM(|A-B|^P).^(1/P)
%   'CITY-BLOCK'   | 'C':  SUM(|A-B|)
%   'COSINE'       | 'O':  1 - (A*B')/||A||*||B||
%   'FOURIER'      | 'F'
%   'TANH'         | 'T':  TANH(P*(A*B')/LENGTH(A) + P);
%   'ANOVA'        | 'A':  ANOVA MODEL
%   'BSPLINE'      | 'B':  BSPLINE MODEL, ORDER P
%   'ANOVABSPLINE' | 'AB': ANOVA-BSPLINE MODEL, ORDER P
%   'ANOVASPLINE1' | 'AS1':ANOVA-SPLINE MODEL, ORDER 1
%   'ANOVASPLINE2' | 'AS2':ANOVA-SPLINE MODEL, ORDER 2
%   'ANOVASPLINE3' | 'AS3':ANOVA-SPLINE MODEL, ORDER 3
%
% In the polynomial case for a non-integer P, the proximity is computed
% as D = SIGN(S+1).*ABS(S+1).^P, in order to avoid problems with negative
% inner products S = A*B'. The features of the objects in A and B may be
% weighted by the weights in the vector WEIGHTS.
%
% SEE ALSO
% PROXM, MAPPINGS, DATASETS

% Copyright: Robert P.W. Duin, r.p.w.duin@prtools.org, Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester



function W = proxxm(A,type,s,weights)

prtrace(mfilename);

if (nargin < 4)
  weights = [];
  prwarning(4,'Weights are not supplied, assuming all 1.');
end

if (nargin < 3) | (isempty(s))
  s = [1 1];
  prwarning(4,'The proximity parameter is not supplied, assuming 1.');
end
if length(s) < 2, s(2) = 1; end

if (nargin < 2) | (isempty(type))
  type = 'd';
  prwarning(4,'Type of the proximity is not supplied, assuming ''DISTANCE''.');
end


% No data, return an UNTRAINED mapping
if (nargin < 1) | (isempty(A))
  W = mapping(mfilename,{type,s,weights});
  W = setname(W,'Proximity mapping');
  return;
end


[m,k] = size(A);

if (isstr(type))
  % Definition of the mapping, just store parameters.
  all = char('polynomial','p','homogeneous','h','exp','e','exp-sum','es','rbf','r',...
             'rbf-sum','rs','sigmoid','s','distance','d','thin-plate','tp',...
             'multiquadric','mq','dsigmoid','ds','minkowski','m','city-block','c',...
             'cosine','o','uniform','u','anova','a','anovaspline1','as1',...
             'anovaspline2','as2','anovaspline3','as3','anovabspline',...
             'ab','bspline','b','spline','fourier','f','tanh','t',...
             'n-thin-plate','ntp');
  if (~any(strcmp(cellstr(all),type)))
    error(['Unknown proximity type: ' type])
  end

  A     = cdats(A,1);
  [m,k] = size(A);
  if isdataset(A)
    W = mapping(mfilename,'trained',{A,type,s,weights},getlab(A),getfeatsize(A),getobjsize(A));
  else
    W = mapping(mfilename,'trained',{A,type,s,weights},[],k,m);
  end
  W = setname(W,'Proximity mapping');


elseif ismapping(type)

  % Execution of the mapping: compute a proximity matrix
  % between the input data A and B; output stored in W.

  w = type;
  [B,type,s,weights] = deal(w.data{1},w.data{2},w.data{3},w.data{4});
  [kk,n] = size(w);
  if (k ~= kk)
    error('Mismatch in sizes between the mapping and its data stored internally.');
  end
  if (~isempty(weights))
    if (length(weights) ~= k),
      error('Weight vector has a wrong length.');
    end
    A = A.*(ones(m,1)*weights(:)');
    B = B.*(ones(n,1)*weights(:)');
  end
  AA = A;
  A  = +A;
  B  = +B;

  switch type
    case {'polynomial','p'}
      D = A*B' + ones(m,n);
      if (s ~= round(s(1)))
        D = sign(D).*abs(D).^s(1);
      elseif (s(1) ~= 1)
        D = D.^s(1);
      end

    case {'homogeneous','h'}
      D = (A*B');
      if (s ~= round(s(1)))
        D = sign(D).*abs(D).^s(1);
      elseif (s(1) ~= 1)
        D = D.^s(1);
      end

    case {'sigmoid','s'}
      D = (A*B');
      D = sigm(D/s(1));

    case {'dsigmoid','ds'}
      D = dist2(B,A).^s(1);
      D = sigm(D/s(2));

    case {'city-block','c'}
      D = zeros(m,n);
      for j=1:n
        D(:,j) = sum(abs(A - repmat(+B(j,:),m,1)),2);
      end

    case {'minkowski','m'}
      D = zeros(m,n);
      for j=1:n
        D(:,j) = sum(abs(A - repmat(+B(j,:),m,1)).^s(1),2).^(1/s(1));
      end

    case {'exp','e'}
      D = dist2(B,A);
      D = exp(-sqrt(D)/s(1));

    case {'exp-sum','es'}
      d = sqrt(dist2(B,A));
      D = sign(s(1))*exp(-d/s(1));
      for z=2:length(s)
        D = D + sign(s(z))*exp(-d/s(z));
      end

      case {'rbf','r'}
        D = dist2(B,A);
        D = exp(-D/(s(1)*s(1)));

      case {'rbf-sum','rs'}
        d = dist2(B,A);
        D = sign(s(1))*exp(-d/(s(1)^2));;
        for z=2:length(s)
          D = D + sign(s(z))*exp(-d/(s(z)^2));
        end

      case {'distance','d'}
        D = dist2(B,A);
        if s(1) ~= 2
          D = sqrt(D).^s(1);
        end

      case {'multiquadric','mq'}
        D = sqrt(dist2(B,A)/s(1) +s(2).^2);

      case {'thin-plate','tp'}
        d = dist2(B,A).^s(1)/s(2);
        D = d.*log(1+d);

      case {'n-thin-plate','ntp'}
        d = dist2(B,A).^s(1)/s(2);
        D = -d.*log(1+d);

      case {'cosine','o'}
        D  = (A*B');
        lA = sqrt(sum(A.*A,2));
        lB = sqrt(sum(B.*B,2));
        D  = 1 - D./(lA*lB');

      case {'tanh','t'}
        D = tanh(s(1)*(A*B')/k + s(1));

      case {'fourier','f'}
        zz = sin(s(1) + 1/2)*2*ones(k,1);
        D  = zeros(m,n);
        for oo=1:m
          for o=1:n
            z = zz;
            I = find(A(oo,:)-B(o,:));
            z(I) = sin(s(1) + 1/2)*(A(oo,I)-B(o,I))./sin((A(oo,I)-B(o,I))/2);
            D(oo,o) = prod(z);
          end
        end

      case {'spline','anovaspline1','as1'}
        D = zeros(m,n);
        for oo=1:m
          for o=1:n
            mAB = min(A(oo,:),B(o,:));
            z   = 1 + A(oo,:).*B(o,:).*(1 + mAB) - ((A(oo,:)+B(o,:))/2).*mAB.^2 + mAB.^3/3;
            D(oo,o) = prod(z);
          end
        end

      case {'anova','a'}
        D = zeros(m,n);
        for oo=1:m
          for o=1:n
            mAB = min(A(oo,:),B(o,:));
            z   = 1 + A(oo,:).*B(o,:).*(1 + 0.5*mAB) - mAB.^3/6;
            D(oo,o) = prod(z);
          end
        end

      case {'bspline','b'}
        D = zeros(m,n);
        for oo=1:m
          for o=1:n
            z = 0;
            AmB = A(oo,:)- B(o,:);
            for r = 0:2*(s(1)+1)
              z = z + (-1)^r*binomial(2*(s(1)+1),r)*(max(0,AmB + s(1)+1 - r)).^(2*s(1) + 1);
            end
            D(oo,o) = prod(z);
          end
        end

      case {'anovaspline2','as2'}
        D = zeros(m,n);
        for oo=1:m
          for o=1:n
            mAB = min(A(oo,:),B(o,:));
            AmB = A(oo,:) - B(o,:);
            AB  = A(oo,:) .* B(o,:);
            ApB = A(oo,:) + B(o,:);
            z   = 1 + AB + AB.^2.*(1 + mAB) - AB.*(ApB).*mAB.^2 + ...
                  (ApB.^2+2*AB).*mAB.^3/3 - 0.5*ApB.*mAB.^4 + 0.2*mAB.^5;
            D(oo,o) = prod(z);
          end
        end

      case {'anovaspline3','as3'}
        D = zeros(m,n);
        for oo=1:m
          for o=1:n
            mAB = min(A(oo,:),B(o,:));
            AmB = A(oo,:)- B(o,:);
            AB  = A(oo,:).*B(o,:);
            ApB = A(oo,:)+B(o,:);
            z   = 1 + AB + AB.^2 + AB.^3.*(1 + mAB) - 1.5 * AB.^2 .*(ApB).*mAB.^2 + ...
                  AB.*(ApB.^2+AB).*mAB.^3 - 0.25 * (3*ApB.^3- 2*B(oo,:).^2 - ...
                  2*A(o,:).^2).*mAB.^4 + 0.6 * (ApB.^2+AB).*mAB.^5 - 0.5 * ApB .*mAB.^6+mAB.^7/7;
            D(oo,o) = prod(z);
          end
        end

      case {'anovabspline','ab'}
        D = zeros(m,n);
        for oo=1:m
          for o=1:n
            z = 0;
            AmB = A(oo,:)- B(o,:);
            for r = 0:2*(s(1)+1)
              z = z + (-1)^r*binomial(2*(s(1)+1),r)*(max(0,BmA + s(1)+1 - r)).^(2*s(1) + 1);
            end
          D(oo,o) = prod(1+z);
        end
      end

      otherwise
        error('Unknown proximity type')
    end
    W = setdat(AA,D,w);

  else
      error('Illegal TYPE argument.')
  end
return;



function D = dist2(B,A)
% Computes square Euclidean distance, avoiding large matrices for a high
% dimensionality data
    m = size(A,1);
    n = size(B,1);
    D = ones(m,1)*sum(B'.*B',1);
    D = D + sum(A'.*A',1)'*ones(1,n);
    D = D - 2 .* (+A)*(+B)';
    % Clean result.
    J = find(D<0);
    D(J) = zeros(size(J));
return
