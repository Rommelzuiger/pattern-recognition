function a = readchicken(norm,cost,dirname)
if nargin < 3
    dirname = ['/data/pr/projects/dissim/data/chickenpieces/chickenpieces_norm' num2str(norm) '.0'];
		if ispc
			dirname = ['s:\root' dirname];
		end
end
if nargin < 2
    dir(dirname)
    a = [];
    return
end
name = ['chickenpieces_norm' num2str(norm) '.0_AngleCostFunction' num2str(cost) '.0.dm'];
a = prtbunke(fullfile(dirname,name));
if length(a)==1 & a==-1
    return
end
a = setprior(a,getprior(a,0));
%a = min(a,a');      % starting from exp4
a = (a + a')/2; 
