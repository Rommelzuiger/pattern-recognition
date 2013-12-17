%SETNAME Mapping for easy name setting
%
%   A = A*SETNAME([],NAME)
%   W = W*SETNAME([],NAME)
%
%Set name of dataset A or mapping W

function a = setname(a,varargin)

if nargin < 1 | isempty(a)
	a = prmapping(mfilename,'combiner',varargin);
else
  if isa(a,'prdataset') || isa(a,'prmapping')
    a = setname(a,varargin);
  else
    error('Illegal input')
  end
end