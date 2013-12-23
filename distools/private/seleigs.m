%SELEIGS Select eigenvalues from a list
%
%   J = SELEIGS(L,ALF)
%
% INPUT
%   L   List of eigenvalues
%   ALF Parameter determining the dimensionality and the eigenvalue-based mapping
%       (0,1)   - fraction of the total (absolute value) preserved variance
%        Inf    - no dimensionality reduction, keeping all dimensions (it's noisy)
%       'p'     - projection into a Euclidean space based on positive eigenvalues only
%       'PARp'  - projection into a Euclidean space based on the PAR fraction of
%                 positive eigenvalues; e.g. ALF = '0.9p'
%       'n'     - projection into a Euclidean space based on negative eigenvalues only
%       'PARn'  - projection into a (negative) Euclidean space based on the PAR fraction
%                 of negative eigenvalues; e.g. ALF = '0.7n'
%       'P1pP2n'- projection into a Euclidean space based on the P1 positive eigenvalues
%                 and P2 negative eigenvalues; e.g. ALF = '0.7p0.1n', ALF = '7p2n'
%       1 .. N  - number of dimensions in total
%       [P1 P2] - P1 dimensions or preserved fraction of variance in the positive subspace
%                 and P2 dimensions or preserved fraction of variance in the negative
%                 subspace; e.g. ALF = [5 10], ALF = [0.9 0.1]
%
% OUTPUT
%   J   Index of selected eigenvalues
%
% DESCRIPTION
% This is a low-level routine for PSEM, KPSEM, KPCA and PEPCA. From a list of
% eigenvalues it selects the ones according to ALF.
%
% SEE ALSO
% MAPPINGS, DATASETS, PSEM, KPSEM, KPCA, PEPCA
%

% Elzbieta Pekalska, e.pekalska@ewi.tudelft.nl
% Faculty of Electrical Engineering, Mathematics and Computer Science,
% Delft University of Technology, The Netherlands


function [J,sig] = seleigs(L,alf,prec)

alfstr = lower(alf);

if strcmp(alfstr,'cut'),
  J = cuteigs(L,prec);
else
  % Check whether alf is a string
  if ~isnumeric(alfstr),
    zp = strfind(alfstr,'p');
    zn = strfind(alfstr,'n');
    if isempty(zp), zp = 0; end
    if isempty(zn), zn = 0; end
    ALF = 0;
    if zp > 0,
      if length(zp) > 1,
        error('Wrong ALF.');
      end
      if zn > zp,
        sstr = alfstr(1:zp-1);
      elseif zn < zp,
        sstr = alfstr(zn+1:zp-1);
      else
        error('Wrong ALF.');
      end
      if ~isempty(sstr),
        ALF = str2num(sstr);
      else
        ALF = length(find(L>0));
      end
    end
    if zn > 0,
      if length(zn) > 1,
        error('Wrong ALF.');
      end
      if zn > zp,
        sstr = alfstr(zp+1:zn-1);
      elseif zn < zp,
        sstr = alfstr(1:zn-1);
      else
        error('Wrong ALF.');
      end
      if ~isempty(sstr),
        ALF = [ALF str2num(sstr)];
      else
        ALF = [ALF length(find(L<0))];
      end
    else
      ALF = [ALF 0];
    end
    alf = ALF;
  end


% alf is now the number of eigenvalues
  if max(size(alf)) > 2 | (alf(1) == 0 & alf(2) == 0),
    error('Wrong ALF.')
  elseif max(size(alf)) == 2,
    J = [];
    Z  = find(L>0);
    if alf(1) > 0,
      if ~isempty(Z),
        J = selnumalf(alf(1),L,Z)';
      else
        prwarning(1,'No positive dimensions are selected.');
      end
    end
    Z  = find(L<0);
    if alf(2) > 0,
      if ~isempty(Z),
        J = [J; selnumalf(alf(2),L,Z)'];
      else
        prwarning(1,'No negative dimensions are selected.');
      end
    end
   else
     J = selnumalf(alf,L,1:length(L));
  end
end


if isempty(J),
  error('Wrong choice of ALF; the subspace cannot be found.');
end

I1  = find(L(J) > 0);
I2  = find(L(J) < 0);
sig = [length(I1) length(I2)];

[ll,K1] = sort(-L(J(I1)));
[ll,K2] = sort(L(J(I2)));

% J - index of selected eigenvalues; SORTED:
% first decreasing positive ones folowed by increasing negative ones
JJ = [];
if sig(1) > 0
  JJ = J(I1(K1));
end
if sig(2) > 0
  JJ = [JJ; J(I2(K2))];
end
J = JJ;
return





function J = selnumalf(alf,L,I)

if alf >= 1,
  if alf == inf,
    cs = cumsum(abs(L(I)));
    nn = min(find(cs./cs(end) == 1));
    J  = I(1:nn)';
  elseif alf > length(I),
    error('ALF is invalid.')
  else
    J = I(1:alf)';
  end
elseif alf > 0,                     % alf in (0,1), percentage
  cs = cumsum(abs(L(I)));
  nn = min(find (cs./cs(end) >= alf));
  J  = I(1:nn)';
elseif alf < 0,
  I   = I(length(I):-1:1);
  alf = abs(alf);
  if alf >= 1,
    if alf > length(I),
      error('ALF is invalid.')
    else
      J  = I(1:alf)';
    end
  else                             % |alf| in (0,1), percentage
    cs = cumsum(abs(L(I)));
    nn = min(find (cs./cs(end) >= alf));
    J  = I(1:nn)';
  end
else
  error('Wrong ALF.');
end
