%PRNIST Load NIST training set in a datafile
%
%	A = PRNIST(CLASS,OBJECTS)
%
% Reads a part of the NIST digit data. CLASS has to be a vector with
% the desired class numbers (between 0 and 9). OBJECTS has to be a
% vector defining the desired object numbers (between 1 and 1000).
% Default OBJECTS = [1:200].

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function a = prnist(n,m)
if nargout == 0
	help(mfilename);
end
if nargin < 1 | isempty(n)
  n = [0:9];
end
if nargin < 2
	m = [1:200];
end
if max(n) > 9 | max(m) > 1000 | min(n) < 0 | min(m) < 1
	error('Class numbers or object numbers out of range')
end

datadir = fileparts(which(mfilename));
if exist(fullfile(datadir,'nisttrain_cell')) ~= 7
	error('The nist data directory could not be found')
end

c = length(n);
a = prdatafile(fullfile(datadir,'nisttrain_cell'));
a = seldat(a,n+1);
J = repmat(m(:),1,c) + repmat([0:c-1]*1000,length(m),1);
a = a(J(:),:);
