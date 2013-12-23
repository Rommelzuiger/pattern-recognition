%MCLASSDC Multi-Class Dissimilarity-based Classifier from Two-Class Discriminants
%
% 	W = MCLASSDC(D,CLASSF,MODE)
%
% INPUT
%   D       Dissimilarity dataset
%   CLASSF  Untrained dissimilarity-based classifier
%   MODE	Way of handling multi-class problems: 'SINGLE' or 'MULTI' 
% 			(optional; default: 'SINGLE')
%
% OUTPUT
%   W       Combined classifier
%
% DEFAULT
%   MODE = 'SINGLE'
% 
% DESCRIPTION
% For the default MODE = 'SINGLE', the untrained classifier CLASSF is called to
% compute C classifiers between each of the C classes in the dataset D and
% the remaining C-1 classes. The result is stored into the combined
% classifier W. There is no combining rule added. The default rule, MAXC
% might be replaced by adding one, e.g. W = W*MEANC.
%
% For the MODE = 'MULTI', the untrained classifier CLASSF is trained between 
% all pairs of classes as well as between each class and all other classes.
% This total set of C^2 classifiers is combined by MINC.	The use of soft labels
% is supported.
% The pairwise classifiers are trained by using a representation set based
% on the traing objects of the corresponding classes only (just in case of
% crisp labels). This is the only difference with the use of MCLASSC. 
%
% EXAMPLE
% W = MCLASSDC(DISTM(GENDATM(100)),QDC([],0.01),'MULTI')
%
% SEE ALSO
% DATASETS, MAPPINGS, MAXC, MINC

% Copyright: R.P.W. Duin, r.duin@ieee.org
% Faculty EWI, Delft University of Technology 


function varargout = mclassdc(D,classf,mode)

prtrace(mfilename);
if nargin < 3, 
	mode = 'single'; 
end
if nargin < 2, 
	classf = []; 
end
if nargin < 1 | isempty(D)
	w = mapping(mfilename,{classf,mode});
	return
end
	
if ~ismapping(classf) | ~isuntrained(classf)
	error('Second parameter should be an untrained mapping.')
end

islabtype(D,'crisp','soft');
isvaldset(D,1,2); 	% At least 1 object per class and two classes

[m,k,c] = getsize(D);
	
varout = {};
if c == 2
	varargout = map(D,classf);
	return
end
lablist = getlablist(D);


switch lower(mode)
	case 'single'
		w = [];
	  %	lablist = getlablist(D);
	  for i=1:c
			if islabtype(D,'crisp')
			  mlab = 2 - (getnlab(D) == i);
			  DD = setlabels(D,mlab);
		  elseif islabtype(D,'soft')
			  Dtargets = gettargets(D);
			  targets  = [Dtargets(:,i) 1-Dtargets(:,i)]; 
			  DD			 = dataset(+D,mlab,targets,'lablist',[1 2]');
		  end
			if nargout(classf.mapping_file) > 1
		 	 	[v,varo] = map(DD,classf);
				varout   = [varout {varo}];
			else
				v = map(DD,classf);
			end
			w = [w,setlabels(v(:,1),lablist(i,:))];
	  end


	case 'multi'
		w = [];
	  nlab = getnlab(D);
		[nl,clab,list] = renumlab(getlablist(D),getfeatlab(D));
	  for z=1:c
		  lab = lablist(z,:);
		  
		  J1 = find(nlab==z);
		  L1 = find(clab==z);
		  if islabtype(D,'crisp')
			  mlab     = ones(m,1);
			  mlab(J1) = zeros(length(J1),1);
			  DD       = setlabels(D,mlab);
		  else
			  problab  = gettargets(D);
			  mlab     = [problab(:,i1) sum(problab,2)-problab(:,z)];
			  DD       = settargets(D,mlab,[1 2]');
		  end		
		  I1 = setdiff([1:c],z);
			if nargout(classf.mapping_file) > 1
		 	 	[v,varo] = map(DD,classf);
				varout   = [varout {varo}];
			else
				v = map(DD,classf);
			end
			w = [w,setlabels(v(:,1),lab)];

		  for t = I1 
			  if islabtype(D,'crisp')
				  J2 = find(nlab==t);
				  L2 = find(clab==t);
				  v  = featsel(k,[L1;L2])*(DD([J1;J2],[L1;L2])*classf);
					%disp([z,t]);
					%parsc(v);
			  else
				  mlab2 = problab(:,[z t]);
				  v     = setlabels(DD,mlab2)*classf;
			  end
			  w = [w,setlabels(v(:,1),lab)];
		  end
	  end
	  w = minc(w);
		
	otherwise
	  error('Unknown mode')
end

w = setname(w,getname(classf));
w = setsize(w,[k,c]);
w = setcost(w,D);
varargout = {w varout};
return
