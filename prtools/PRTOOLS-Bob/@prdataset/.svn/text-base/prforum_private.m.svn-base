

function prforum_private(a,k)

if isstr(k), k = str2num(k); end

k = floor(k);

if k >= 1000000 | k < 1
	n = 0;
else
	if exist('prtools/@dataset/dataset') == 0
		error('PRTools not found');
	end
	seed = rand('state');
	rand('state',893267);
	L = randperm(999999);
	rand('state',seed);
	n = L(k);
end
disp(n)
