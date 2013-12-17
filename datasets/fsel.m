% FSEL Feature selection
%   [W,LIST] = FSEL (DATA, ALGORITHM, CRITERION, P) selects P features
%   out of the D in dataset DATA which give the highest CRITERION value. The
%   search algorithm ALGORITHM can be either 'individual' (best-P), 
%   'forward', 'backward', '+l-r' (plus-l-takeway-r) or 'b&b' (branch & bound). 
%   It returns a mapping W and a feature LIST ordered by criterion.
%
%   If P is not given, all D features are ranked. 
%   For possible values for CRITERION, see feateval. 
%   Defaults: ALGORITHM = 'individual', CRITERION = 'NN'.

function [W,list] = fsel (data, algorithm, criterion, k, l, r)

	[n,m] = size(data);
	
	if (nargin < 6), r = 0; 														end;
	if (nargin < 5), l = 1; 														end;
	if (nargin < 4), k = m;       											end;
	if (nargin < 3), criterion = 'NN'; 									end;
	if (nargin < 2), algorithm = 'individual';					end;
	if (nargin < 1), error ('insufficient arguments'); 	end;

	switch algorithm
		case {'i', 'individual'},
			[W,R] = featseli (data, criterion, k);
			list = R(1:k,3);
		case {'f', 'forward'},
			[W,R] = featself (data, criterion, k);
			list = R(1:k,3);
		case {'b', 'backward'},
			[W,R] = featselb (data, criterion, 1);
			R(:,3) = -R(:,3);
			last_feature = setdiff((1:m)',R(:,3));
			R = [R; 0 -1 last_feature];
			list = R(end:-1:end-k+1,3);
			[W,R] = featselb (data, criterion, k);
		case {'l', '+l-r' },
			[W,R] = featsellr (data, criterion, k, l, r);
			list = (+W)';
		case {'bnb', 'b&b' },
			W = featselo (data, criterion, k);
			list = (+W)';
%		case {'o', 'optimal'},
%			j = 1; prev_features = 1:m;
%			for i = m-1:-1:1
%				W = featselo(data, criterion, i);
%				dropped_feature = setdiff (prev_features, +W);
%				R(j,:) = [ j -1 -dropped_feature ];
%				j = j + 1;
%				prev_features = +W;
%			end;
%			R(:,3) = -R(:,3);
%			list = R(end:-1:end-k,:);
		otherwise,
			error ('unknown algorithm');
	end;

%	W = mapping('featsel',list,list,m,k,1,[]);
	W = featsel(m,list);

return
