%NIST_EVAL Evaluation of NIST digits classifier
%
%		E = NIST_EVAL(DIG2DATA,W,N)
%
% INPUT
%   DIG2DATA  M-File to convert digits images into a dataset
%   W         Classifier
%   N         Number of digits per class to be used (default 100);
%
% OUTPUT
%   E         Estimate of classification error
%
% DESCRIPTION
%   This routine reads N digits per class (0:9) into a datafile.
%   DIG2DATA is used to convert them into a dataset.  DIG2DATA should
%   be given by the name (string) of an m-file in the search path
%   (without the .m extension!). The resulting dataset is given to the
%   classifier W to test it.

function e = nist_eval(file,w,n)

if exist(file)  ~= 2
	error(['File ' file ' not found in search path'])
end

ismapping(w);

if nargin < 3
	n = 100;
	warning('Assuming 10 test objects per class')
end

if n < 1 | n > 100
	error('Number of test objects per class should be between 1 and 100')
end

datdir=which('nist_eval');
datdir=fileparts(datdir);
s = prload(fullfile(datdir,'.in4085','.in4085'),'-mat');
m = s.b;
m = set(m,'rootpath',fullfile(datdir,'.in4085'));
%m = datafile(fullfile(datdir,'.in4085'),'cell');
%m = gendat(m,ones(1,10));
R = randperm(10*n);
m = m(R,:);
nlab = getnlab(m);
m = setnlab(m,zeros(10*n,1));
a = feval(file,m);
if ~isdataset(a)
	error(['Command ' file ' does not return a proper dataset']);
end
a = setlablist(a,getlablist(m));
a = setnlab(a,nlab);
a = setprior(a,0);
e = a*w*testd;
