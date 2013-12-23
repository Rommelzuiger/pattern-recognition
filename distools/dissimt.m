%DISSIMT  Fixed DISsimilarity-SIMilarity transformation
%
%    W = DISSIMT([],TYPE,PAR)
%           OR
%    T = DISSIMT(P,TYPE,PAR)
%
% INPUT
%   P     NxN or NxM proximity matrix or dataset
%   TYPE  Type of the transformation that acts on proximity data
%         (optional; default: 'POW'):
%         'I',    'identity':   P
%         'SCALE', 'scale': 		P ./ PAR
%         'NEG',  'negative':   -P
%         'NEG2', 'negative2':  -P.^2
%         'REV',  'reverse':    1./(PAR + P)
%         'REV2', 'reverse':    1./(PAR^2 + P.^2)
%         'D2S',  'S2D':        1-P
%         'D2SSQ','S2DSQ':      sqrt(1-P)
%         'POW',  'power':      P.^PAR
%         'EXP',  'exponent':   exp(-P./PAR)
%         'EXP2', 'exponent2':  exp((-P.^2)./(PAR^2))
%         'GAUSS':              exp((-P.^2)./(2*PAR^2))
%         'LOG',  'logarithm':  log(P/PAR)
%         'NLOG', 'neglog':     -log((P+1)./PAR)
%         'NLOG2','neglog2':    -log((P.^2+1)./PAR^2)
%         'SIGM', 'sigmoid':    2./(1+exp(-P./PAR))-1
%         'SIGM2','sigmoid2':   2./(1+exp(-(P.^2)./PAR^2))-1
%         'SIM2DIS'
%         'DIS2SIM'
%   PAR   Parameter, PAR > 0 (optional; default: 1)
%
% OUTPUT
%   W     Fixed mapping
%   T     NxN or NxM proximity matrix or dataset
%
% DESCRIPTION
% A fixed mapping that transforms proximity data. Note that in some cases
% similarity is transformed into dissimilarity, or vice versa.
% If a proximity data P is supplied, the data are directly mapped, so there
% is no mapping defined. Hence, DISSIMT([],TYPE,PAR) is a mapping, while
% DISSIMT(P,TYPE,PAR) returns data.
%
% DEFAULT
%   PAR = 1
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function W = dissimt(D,type,p)
if nargin < 3,
  p = 1;
end
if nargin < 2,
  type = 'pow';
end
if nargin == 0 | isempty(D)
  W = mapping(mfilename,'fixed',{type,p});
  W = setname(W,'Fixed proximity transformation');
  return
end

if p <= 0 or p == inf,
  error('Wrong PAR.');
end

switch lower(type)
case {'i','identity'}
   W = D;
case {'scale'}
   W = D./p;
case {'neg','negative'}
   W = -D;
case {'neg2','negative2'}
   W = -D.^2;
case {'rev','reverse'}
   W = 1./(p+D);
case {'rev2','reverse2'}
   W = 1./(p^2+D.^2);
case {'d2s','s2d'}
   if any(D > 1 | D < 0),
     error('Proximity data should be in [0,1].')
   end
   W = 1-D;
case {'d2ssq','s2dsq'}
   if any(D > 1 | D < 0),
     error('Proximity data should be in [0,1].')
   end
   W = sqrt(1-D);
case {'pow','power'}
   W = D.^p;
case {'exp','exponent'}
   W = exp(-D/p);
case {'exp2','exponent2'}
   W = exp(-D.^2/p^2);
case {'gauss'}
   W = exp(-D.^2/(2*p^2));
case {'log','logarithm'}
   W = log(D/p);
case {'nlog','neglog'}
   W = -log((D+1)/p);
case {'nlog2','neglog2'}
   W = -log((D.^2+1)/p^2);
case {'sigm','sigmoid'}
   W = 2./(1+exp(-D/p))-1;
case {'sigm2','sigmoid2'}
   W = 2./(1+exp(-D.^2/p^2))-1;
case {'dis2sim'}
	m = size(D,1);
	one = ones(m,1);
	H = eye(m) - one*one'/m;
	W = -H*D.^2*H/2;
case {'sim2dis'}
	m = size(D,1);
	DD = diag(+D);
	W = sqrt(-2*D+repmat(DD,1,m)+repmat(DD',m,1));
otherwise
   error('Transformation is unknown.');
end
return
