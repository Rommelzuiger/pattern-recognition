%ISPE_EM Test on pseudo-Euclidean mapping
%
%  N = ISPE_EM(W)
%      ISPE_EM(W)
%
% INPUT
%  W    input mapping
%
% OUTPUT
%  N    logical value
%
% DESCRIPTION
% Returns true for pe_em mappings. If no output is required,
% false outputs are turned into errors. This may be used for
% assertion.
%
% SEE ALSO
% ISPE_EM

function n = ispe_em(w)

	prtrace(mfilename);
	
	if isa(w,'mapping') & strcmp(w.mapping_file,'pe_em')
		n = 1;
	else
		n = 0;
	end

	% generate error if input is not a psem mapping
	% AND no output is requested (assertion)

	if nargout == 0 & n == 0
		error([newline '---- pe_em mapping expected -----'])
	end

return
