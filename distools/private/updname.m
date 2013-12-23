%UPDNAME Update of a distance/similarity dataset
%
%   NAME = UPDNAME (NAME)
%
% A low level routine used in other functions

% Copyright: Elzbieta Pekalska, e.pekalska@ewi.tudelft.nl
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands



function name = updname(name)

str = lower(name);
NAMES = {'distance matrix', 'similarity matrix', 'kernel matrix'};

ok = 0;
z  = 1;
while ~ok & z<=length(NAMES)
  NAME = NAMES{z};
  m = length(NAME);
  k = findstr(str,NAME);
  if ~isempty(k),
    if length(str) >=k+m+3 & strcmp(str(k+m+1:k+m+3),'for'),
      name(k:k+m+3)=[];
    else
      name(k:k+m-1)=[];
    end
    if ~isempty(name) & strcmp(name(1), ' '),
      name(1)= [];
    end
    ok = 1;
  end
  z = z+1;
end
if isempty(name),
  name = 'Data';
end
