%PE_AFFINE Define and execute Pseudo Eclidean Affine mappings
%
%   W = PE_AFFINE(W,SIG)
%   F = PE_AFFINE(A,W)
%   F = A*W

function argout = pe_affine(arg1,arg2)

if ismapping(arg1)  % definition of PE mapping
	isaffine(arg1);
	if length(arg2) ~= 2
		error('Signature should have two components')
	end
	[n1,n2] = size(arg1);
	if sum(arg2) ~= n1
		error('Sum of signature numbers should equal input dimension of mapping')
	end
	par.map = arg1;
	par.sig = arg2;
	argout = mapping(mfilename,'trained',par,getlabels(arg1),n1,n2);
else
	isdataset(arg1);
	ismapping(arg2);
	w = getdata(arg2,'map');
	sig = getdata(arg2,'sig');
	if sig(2) ~= 0;
		data = getdata(w);
		data.rot = -data.rot;
		data.offset = zeros(1,size(arg2,2));
		v = setdata(w,data);
		[m,k] = size(arg1);
		if sig(1) ~= 0
			argout = [arg1(:,1:sig(1)) zeros(m,sig(2))]*w;
		else
			argout = setdata(arg1,zeros(m,k))*w;
		end
		argout = argout + [setdata(arg1,zeros(m,sig(1))) arg1(:,sig(1)+1:end)]*v;
	else % no negative directions: back in Euclidean case
		argout = arg1*w;
	end
end
