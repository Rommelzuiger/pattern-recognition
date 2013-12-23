%DLPC LP-classifier on dissimilarity (proximity) data
%
%   [W1,W2,W3] = DLPC(D,BIAS,TYPE,PARAM)
%
% INPUT
%   D     Dissimilarity (proximity) dataset
%   BIAS  YES or NO (optional; default: 1 (YES))
%   TYPE  Type of a classifier
%         'SIMPLE'   - the most simple formulation; no sparse solution; PARAM = [];
%         'STANDARD' - minimization of the training misclassification errors;
%                      no sparse solution; PARAM = [];
%         'C-SPARSE' - sparse solution; a formulation similar to the LP_1 SVM;
%                       PARAM is a tradeoff parameter, similar as in the traditional
%                       SVM; (optional; DEFAULT: 1).
%         'MU-SPARSE' - sparse solution; a formulation similar to the LP_1 SVM,
%                       based on the paper of Graepel, Herbrich, Smola etc
%                       'Classification on proximity data with LP-machines'.
%                       PARAM is a tradeoff parameter, usually PARAM = 0.05 or 0.1.
%                       It is an upper bound on the misclassfied training objects.
%                       So, for well separable problems, PARAM = 0.01 or PARAM = 0.02.
%                       (optional; DEFAULT: the LOO 1-NN error * 1.3).
%  PARAM Parameter connected to the TYPE, as above
%
% OUTPUT
%   W1  LP-Classifier in the complete dissimilarity space
%   W2  LP-Classifier in a reduced dissimilarity space
%   W3  Object selection mapping; the indices of support objects are in +W3.
%
% DEFAULTS
%   BIAS  = 1
%   TYPE  = 'STANDARD'
%   PARAM = []
%
% DESCRIPTION
% Classification problem on a N x M dissimilarity dataset D with LP-machines.
% D should be described by both label and feature lists. If D is a square,
% symmetric matrix, then the feature list should be the same as the label list.
%
% Assume a 2-class problem. Let DLPC select J support objects. Then:
% W1 is an M x 2 classifier in the original dissimilarity space, W2 is an J x 2
% classifier in the dissimilarity space defined by the J support objects
% and W3 is an M x R feature selection such that W1 = W3 * W2.
% Note that the indices of the support objects can be retrieved by +W3.
%
% A linear classifier is built on D:
%
%     f(D(x,*)) = diag(Y) * D(x,*) * W + W0,
%
% where Y are labels (+/- 1) and W are the weights. If BIAS is 1, then W0 is also
% sought, otherwise it equals 0, hence the hyperplane is forced to go through the origin.
%
% For C-class problems, C classifiers are trained, one against all others.
% In such a case, only W1 is returned and W3 in now NOT a feature selection,
% but directly the indices of the support objects.


%
% Let [n,k] = size(D). Assume a two-class problem.
% Any multi-class problem is converted one-against-all to two-class problems.
% Y are the labels (converted +/-1)
% D_Y = diag(Y_r) * D * diag(Y_c),
% where Y_r are the labels for rows of D and
% Y_c are the labels for columns.
% alpha is the sought solution.
%
% Internal - the optimization schema (A,b,f) refer to the constraints in
% the standard Matlab LINPROG procedure.
%
% 'simple':
%   min    0^T * alpha(1:k) (= 0)
%   s.t.   [D_Y  Y]  * [alpha(1:k) alpha_0)] >= 1
%
%   A = [D_Y  Y]  is  n x (k+1)
%   b = 1(n,1)
%   f = 0(k+1,1)
%
%
% 'standard':
%   min    p^T * beta(1:n),    p = 1/n_i for the class_i
%   s.t.   [D_Y Y] * [alpha(1:k) alpha_0] + beta(1:n) >= 1,
%          beta >= 0
%
%   A = [D_Y  Y  eye(n)]
%   b = 1(n,1)
%   f = [0(k+1,1) p]^T
%
%
% 'c-sparse':
%   min    sum(alpha(1:k)) + sum (alpha2(1:k)) + C*sum(ksi)
%   s.t.   [D_Y Y] * [alpha(1:k)-alpha2(1:k)  alpha_0] + ksi(1:n) >= 1,
%          ksi >= 0
%          alpha, alpha2 >= 0
%
%   A = [D_Y  -D_Y Y  eye(n,n)]
%   b = 1(n,1)
%   f = [0(k+1,1)]^T
%
%
% 'mu-sparse':
%   min    sum(ksi)/n - mu*g
%   s.t.   [D_Y -D_Y Y] * [alpha(1:k)-alpha2(1:k)  alpha_0] + ksi(1:n) >= g,
%          sum(alpha) + sum(alpha2) = 1
%          alpha, alpha2 >= 0
%          ksi >= 0
%
%   A = [D_Y  -D_Y  Y  eye(n,n)]
%   b = 1(n,1)
%   f = [0(2k+1,1)]
%


% Elzbieta Pekalska, Robert P.W. Duin, e.pekalska@tudelft.nl
% Faculty of Electrical Engineering, Mathematics and Computer Science,
% Delft University of Technology, The Netherlands.



function [W1,W2,W3] = dlpc (d,is_w0,type,par,usematlab,prec)
if nargin < 6, prec = 1e-7; end
if nargin < 5, usematlab = 0; end
if nargin < 3 | isempty(type), type =  'standard'; end
if nargin < 4
    par = [];
end
if nargin < 2 | isempty(is_w0), is_w0 = 1; end
if nargin < 1 | isempty(d)
  W1 = mapping(mfilename,{is_w0,type,par,usematlab});
  W1 = setname(W1,'DLPC');
  W2 = [];
  W3 = [];
  return
end


if ~isdataset(d),
  error('The first parameter should be a dataset.')
end
if ~isnumeric(is_w0) | (is_w0 ~= 0 & is_w0 ~= 1),
  error('The second parameter should be 0 or 1.');
end

if isempty(par),
  switch upper(type)
    case 'MU-SPARSE',
      par = max(0.01,1.3*testkd(d,1,'loo')); % upperbound error: 1.3 * loo 1-nn error
    case 'C-SPARSE',
      par = 1;
    case {'SIMPLE','STANDARD'},
      par = [];
    otherwise
      error('Wrong type.')
  end
end

lab      = getnlab(d);
lablist  = getlablist(d);
featlab  = getfeatlab(d);
[m,k,C]  = getsize(d);

[nl, fe, fl] = renumlab(lablist,featlab);
if max(fe) > size(lablist,1)
  error('Feature labels of the dataset do not match with class labels.')
end



z = (is_w0 > 0);  % is bias used or not?

% This is the status of the optimization procedure.
% For GLPK, this is the exit code; see GLPKMEX for details.
% For Matlab LINPROG, if negative then no solution is found.

status = 1;


% If more than 2 classes, do one against all others.
if C > 2,
% W1 = mclassc(d,mapping(mfilename,{is_w0,type,par,usematlab}));

  W1 = [];
  W2 = [];
  W3 = [];
  N = [];
  for i=1:C
    mlab = 2 - (lab == i);
    mfe  = 2 - (fe == i);
    dd   = setlabels(d,mlab);
    dd   = setfeatlab(dd,mfe);
    if ~isempty(d.prior)
      dd = setprior(dd,[d.prior(i),1-d.prior(i)]');
    end
    [v1,v2,v3]= dlpc(dd,is_w0,type,par,usematlab);
    j = +v3;
    if isempty(v1),
      W1 = [];
      W2 = [];
      W3 = [];
      prwarning(1,'No solution found.');
      return;
    end
    W1 = [W1,setlabels(v1(:,1),lablist(i,:))];
    W2 = [W2;setlabels(v2(:,1),lablist(i,:))];
    W3(j) = ones(length(j),1);
    N = [N j];
  end
  [N1,N2,N3] = unique(N);
  W3 = featsel(k,N1);
  %disp(size(W3,2))
  W2 = featsel(length(N1),N3)*W2;
  return

else

  Y1 = 3 - 2 * lab;   % labels     +/-1
  Y  = 3 - 2 * fe;    % featlabels +/-1

  alpha(1:k+1,1) = 0;


  switch type
    case 'simple',
      f = zeros(k+z,1);
      b = -ones(m,1);
      if is_w0,
        A = -[(Y1*Y').* +d  Y1];
      else
        A = -[(Y1*Y').* +d];
      end

%     if (exist('glpkmex')>0) & (usematlab==0)
%       smin     = 1;  % solve minimum
%       ctype    = char(ones(m,1)*abs('U'));    % Sign of inequalities
%       vartype  = char(ones(k+z,1)*abs('C'));  % Continous variables
%       lpsolver = 2;                           % Interior Point Method
%       lpsolver = 1;                           % Revised Simlex Method
%       params.msglev = 0; % no outputs
%       [al,fval,status] = glpkmex(smin,f,A,b,ctype,[],[],vartype,params,lpsolver);
%     else
%       [al,fval,status] = linprog(f,A,b);
%     end

      [al,fval,status] = linprog(f,A,b);
      alpha(1:k+z) = al;



    case 'standard',
      L    = ones(k,1);
      I    = find(Y==1);
      if ~isempty(I)
        L(I) = L(I)/length(I);
      end
      J    = find(Y==-1);
      if ~isempty(J)
        L(J) = L(J)/length(J);
      end

      f  = [zeros(k+z,1); L];
      lb = [-Inf .*ones(k+z,1); zeros(k,1)];
      ub = Inf .* ones(2*k+z,1);
      b  = -ones(m,1);
      if is_w0,
        A  = -[(Y1*Y').* +d  Y1  eye(m,k)];
      else
        A  = -[(Y1*Y').* +d  eye(m,k)];
      end

%     if (exist('glpkmex')>0) & (usematlab==0)
%       smin     = 1;  % solve minimum
%       ctype    = char(ones(m,1)*abs('U'));     % Sign of inequalities
%       vartype  = char(ones(2*k+z,1)*abs('C')); % Continous variables
%%        lpsolver = 2;                          % Interior Point Method
%       lpsolver = 1;                            % Revised Simlex Method
%       params.msglev = 0; % no outputs
%       [sss,hostname] = unix('hostname');
%       hostname = hostname(1:end-1);
%       if strcmp(hostname,'saturnus') | strcmp(hostname,'polaris') | strcmp(hostname,'neptunus')
%         [al,fval,status] = glpkmex_redhat(smin,f,A,b,ctype,lb,ub,vartype,params,lpsolver);
%       else
%         [al,fval,status] = glpkmex(smin,f,A,b,ctype,lb,ub,vartype,params,lpsolver);
%       end
%     else
%       [al,fval,ststus] = linprog(f,A,b,[],[],lb,ub);
%     end

    [al,fval,ststus] = linprog(f,A,b,[],[],lb,ub);
      alpha(1:k+z) = al(1:k+z);



    case 'c-sparse',
      L  = ones(k,1);
      ub = Inf .* ones(3*k+z,1);
      lb = [zeros(2*k,1); -Inf.*ones(z,1); zeros(k,1)];
      b  = -ones(m,1);

      dd = (Y1*Y').* +d;
      if is_w0,
        f  = [ones(2*k,1); 0; par*L];
        A  = -[dd  -dd  Y1 eye(m,k)];
      else
        f  = [ones(2*k,1); par*L];
        A  = -[dd  -dd  eye(m,k)];
      end
      if (exist('glpkmex')>0) & (usematlab==0)
        smin     = 1;  % solve minimum
        ctype    = char([ones(m,1)*abs('U')]);       % Sign of inequalities
        vartype  = char(ones(3*k+z,1)*abs('C'))  ;   % Continuous variables
%         lpsolver = 1;                               % Revised Simlex Method
        lpsolver = 2;                               % Interior Point Method
        params.msglev = 0;    % no outputs
        params.itlim  = 400;  % iteration limit
        [sss,hostname] = unix('hostname');
        hostname = hostname(1:end-1);
        if strcmp(hostname,'saturnus') | strcmp(hostname,'polaris') | strcmp(hostname,'neptunus')
          [al,fval,status] = glpkmex_redhat(smin,f,A,b,ctype,lb,ub,vartype,params,lpsolver);
        else
          [al,fval,status] = glpkmex(smin,f,A,b,ctype,lb,ub,vartype,params,lpsolver);
        end
      else
        [al,fval,status] = linprog (f,A,b,[],[],lb,ub);
      end
      alpha(1:k) = al(1:k) - al(k+1:2*k);
      if is_w0,
        alpha(k+1) = al(2*k+1);
      end



    case 'mu-sparse',
      L   = ones(k,1)/k;
      f   = [zeros(2*k+z,1); L; -par];
      ub  = Inf .* ones(3*k+1+z,1);
      lb  = [zeros(2*k,1); -Inf.*ones(z,1); zeros(k+1,1)];
      Aeq = [ones(2*k,1); zeros(k+1+z,1)]';
      beq = 1;
      b   = zeros(m,1);
      dd  = (Y1*Y').* +d;

      if is_w0,
        A  = -[dd  -dd  Y1  eye(m,k) -ones(m,1)];
      else
        A  = -[dd -dd  eye(m,k) -ones(m,1)];
      end

      if (exist('glpkmex')>0) & (usematlab==0)
        smin     = 1;  % solve minimum
        ctype    = char([ones(m,1)*abs('U'); 'S']);  % Sign of inequalities
        vartype  = char(ones(3*k+1+z,1)*abs('C'));   % Continous variables
%        lpsolver = 1;                                % Revised Simlex Method
        lpsolver = 2;                               % Interior Point Method
        params.msglev = 0;    % no outputs, but doesn't seem to work
        params.itlim  = 400;  % iteration limit
        [sss,hostname] = unix('hostname');
        hostname = hostname(1:end-1);
        if strcmp(hostname,'saturnus') | strcmp(hostname,'polaris') | strcmp(hostname,'neptunus')
          [al,fval,status] = glpkmex_redhat(smin,f,[A; Aeq],[b; beq],ctype,lb,ub,vartype,params,lpsolver);
        else
          [al,fval,status] = glpkmex(smin,f,[A; Aeq],[b; beq],ctype,lb,ub,vartype,params,lpsolver);
        end
      else
        [al,fval,status] = linprog(f,A,b,Aeq,beq,lb,ub);
      end

      alpha(1:k) = al(1:k) - al(k+1:2*k);
      if is_w0,
        alpha(k+1) = al(2*k+1);
      end

    otherwise
      error ('Wrong type.');
    end
  end

  if (status <= 0) | (status > 181 | status == 150),
    prwarning(1,'Fisher classifier is trained.');
    W1 = fisherc(d);
    W2 = W1;
    W3 = featsel(k,[1:k]);
    return;
  end


  % Choose support objects
  ss = sum(abs(alpha(1:k)));
  J  = find(abs(alpha(1:k)) > ss*prec);
  if ~isempty(J),
    W3 = featsel(k,J);
    w = [Y; 1] .* alpha(1:k+1);
    W2 = affine(w(J),w(k+1),d(:,J),lablist,k,2);
    W2 = cnormc(W2,d(:,J));
    W1 = W3*W2;
    W1 = setname(W1,'DLPC');
    W2 = setname(W2,'DLPC');
  else
    prwarning(1,'No support objects found. Fisher classifier is trained.');
    W1 = fisherc(d);
    W2 = W1;
    W3 = featsel(k,[1:k]);
    return;
  end
  % disp(size(W3,2))
return;
