%DISNORM Normalization of a dissimilarity matrix
%
% 	[V,DN] = disnorm(D,OPT)
%
% INPUT
%   D 	NxN or NxM dissimilarity matrix or dataset
% 	OPT 'max':  maximum dissimilarity is set to 1 by global rescaling
% 		'mean': average dissimilarity is set to 1 by global rescaling (default)
%
% OUTPUT
%   V 	Fixed mapping
%   DN	Normalized dissimilarity data
%
% DEFAULT
%   OPT = 'mean'
%

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester


function [V,D] = disnorm(D,opt)
if nargin < 2,
	opt = 'mean';
end

if nargin == 0 | isempty(D)
  V = mapping(mfilename,'fixed',{opt});
  V = setname(V,'DISNORM');
  return
end

%DEFINE mapping
if isstr(opt)
	opt = lower(opt);
	if strcmp(opt,'mean')
		m = mean(setdiff(D(:),D(1:size(D,1)+1:end)));
	elseif strcmp(opt,'max')
		m = max(D(:));
	else
		error('Wrong OPT.')
	end
	if nargout > 1
		D = D./m;
	end	
    V = mapping(mfilename,'fixed',{m},[],size(D,2),size(D,2));
	return;
end


% APPLY mapping
if isnumeric(opt)
	V = D./opt;
	return;
end			
	
