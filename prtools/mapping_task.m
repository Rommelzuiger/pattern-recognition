%MAPPING_TASK Find task of actual mapping call
%
%   STATE = MAPPING_TASK(ARGIN,TASK)
%
% INPUT
%   ARGIN   - Arguments of a mapping call
%   TASK    - 'definition','fixed execution','training',
%             'trained execution'

% OUTPUT
%   STATE   - TRUE or FALSE
%
% DESCRIPTION
% This routine facilitates the parsing of calls to PRTools mapping
% routines.
%
% SEE ALSO 
% DATASETS, MAPPINGS, MAPPING_TASK, TRAINED_CLASSIFIER, TRAINED_MAPPING

% Copyright: Robert P.W. Duin, prtools@rduin.nl

function state = mapping_task(argin,task)

if nargin < 2
  error('No task specified')
end

state = false;
switch lower(task)
  case {'definition'}
    if isempty(argin) || isempty(argin{1})
      state = true;
    end
  case {'fixed execution','fixed_execution','fixedexecution'}
    if (nargin > 0) && (isdouble(argin{1}) || isa(argin{1},'prdataset')) && ...
        (numel(argin) == 1 || ~ismapping(argin{2}))
      state = true;
    end
  case ('training')
    if (nargin > 0) && (isdouble(argin{1}) || isa(argin{1},'prdataset')) && ...
        (numel(argin)==1 || ~ismapping(argin{2}) || isuntrained(argin{2}))
      state = true;
    end
  case {'execution','trained execution','trained_execution','trainedexecution'}
    if (nargin > 1) && (isdouble(argin{1}) || isa(argin{1},'prdataset')) && ...
        ismapping(argin{2})
      state = true;
    end
  otherwise
    error('Wrongly specified mapping task')
end

return
    