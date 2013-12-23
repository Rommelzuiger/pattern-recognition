%CHICKENPIECES
%
%   D = CHICKENPIECES(NORM,COST,DIR)
%
% Reads the original, asymmetric distances

function a = chickenpieces(norm,cost,dirname)
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
%a = (a + a')/2; 
a = setname(a,['Chickenpieces-' int2str(norm) '-' int2str(cost)]);
desc = ['This is one of the chickenpieces dissimilarity matrices as made available by Bunke et.al. ' ...
  'Every entry is a weighted edit distance between two strings representing the contours of 2D blobs. ' ...
  'Contours are approximated by vectors of length ' num2str(norm) ', using their angles are the replacement costs. ' ...
  'The costs for insertion and deletion are ' num2str(cost)];
ref = ['H. Bunke, H., U. Buhler, Applications of approximate string matching to 2D shape recognition, ' ...
  'Pattern recognition 26 (1993) 1797-1812'];
link = {{'The original datasets', 'http://www.iam.unibe.ch/fki/databases/string-edit-distance-matrices'}, ...
  {'The PRTools versions','http://prtools.org/files/chickenpieces.zip'}};
a = setuser(a,desc,'desc');
a = setuser(a,ref,'ref');
a = setuser(a,link,'link');
