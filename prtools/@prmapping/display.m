%DISPLAY Display mapping information

% $Id: display.m,v 1.4 2009/02/28 17:32:41 duin Exp $

function out=display(w,space)
			if nargin < 2, space = ''; end
	[k,c] = size(w);
	kk = num2str(k);
	cc = num2str(c);

	wname = getname(w);
  if isempty(w)
    name = 'Empty (unity) ';
	elseif ~isempty(wname)
		name = [wname ', '];
	else
		name = '';
	end

	d = w.mapping_file; 
	s = w.mapping_type;
	if k*c == 0
		;
	else
		s = [kk ' to ' cc ' ' s];
	end

	if isempty(d)
		s = [space name s '  mapping'];
	elseif isclassifier(w)
		s = [space name s ' classifier --> ' d];
	else
		s = [space name s '  mapping   --> ' d];
	end

	if nargout == 0
		disp(s);
	else
		out = s;
	end
return
