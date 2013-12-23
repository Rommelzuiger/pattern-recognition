%ISPSEM Test on pseudo-Euclidean mapping
%
%  N = ISPSEM(W)
%      ISPSEM(W)
%
% INPUT
%  W    input mapping
%
% OUTPUT
%  N    logical value
%
% DESCRIPTION
% Returns true for psem mappings. If no output is required,
% false outputs are turned into errors. This may be used for
% assertion.
%
% SEE ALSO
% PSEM

function n = ispsem(w)

	prtrace(mfilename);
	
	if isa(w,'mapping') & strcmp(w.mapping_file,'psem')
		n = 1;
	else
		n = 0;
	end

	% generate error if input is not a psem mapping
	% AND no output is requested (assertion)

	if nargout == 0 & n == 0
		error([newline '---- psem mapping expected -----'])
	end

return
