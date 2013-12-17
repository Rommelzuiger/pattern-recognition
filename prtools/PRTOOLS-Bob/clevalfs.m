%CLEVALFS Classifier evaluation (feature definition, feature size curve)
% 
%   E = CLEVALFS(A,FEATDEF,CLASSF,FEATSIZES,LEARNSIZE,NREPS,S,TESTFUN)
% 
% INPUT
%   A          Training dataset.
%   FEATDEF    Untrained mapping(s), defining the feature space
%              (possibly cell array)
%   CLASSF     Untrained classifier(s) to be tested (possibly cell array)
%   FEATSIZES  Vector of feature sizes (default: all sizes)
%   LEARNSIZE  Number of objects/fraction of training set size
%              (see GENDAT)
%   NREPS      Number of repetitions (default: 1)
%   S          Independent test dataset (optional)
%   TESTFUN    Mapping,evaluation function (default classification error)
%
% OUTPUT
%   E          Cell array [#FEATDEF,#CLASSF] with results (structures)
%              See PLOTE for a description
%
% DESCRIPTION
% This routine is an extension of CLEVALF which makes a feature curve of a
% given dataset. CLEVALFS produces feature curves for every feature
% representation defined in one of the cells of FEATDEF (e.g. based on PCAM
% or FEATSELS). The result is a set of feature curves that can be plotted
% by PLOTE.
%
% The following steps are taken:
% - Split data dataset A in a trainset and and testset S according to
% LEARNSIZE.
% - Compute for an element in FEATDEF a feature representation. This
% includes a ranking of the features.
% - Compute for this set of features and one of the classifiers in CLASSF a
% feature curve for the feature sizes as defined by FEATSIZES.
% - repeat the last line for all classifiers in CLASSF.
% - repeat the last three lines for all feature definitions in FEATDEF.
% - Store all results as cells in E.
% 
% This function uses the RAND random generator and thereby reproduces only
% if its seed is saved and reset.
%
% SEE ALSO 
% MAPPINGS, DATASETS, CLEVAL, CLEVALF, TESTC, PLOTE, GENDAT

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function e = clevalfs(a,featdef,classf,featsizes,learnsize,n,s,testfun)
	
  if (nargin < 8), testfun = []; end
  if (nargin < 7), s = []; end
  if (nargin < 6), n = 1; end
  if (nargin < 5), learnsize = 0.5; end
  if (nargin < 4), featsizes = []; end
  if (nargin < 3)
    error('Insufficient parameters specified')
  end;
  
  if ~iscell(classf), classf = {classf}; end
  if ~iscell(featdef), featdef = {featdef}; end
  isdataset(a);
  a = setprior(a,getprior(a));
  
  if numel(classf) > 1
    e = cell(numel(featdef),1);
  else
    e = cell(1,numel(featdef));
  end
  
  if isempty(s)
    if ismapping(learnsize)
      [t,s] = a*learnsize;
    else
      [t,s] = gendat(a,learnsize);
    end
  else
    isdataset(s);
    t = a;
  end
  
  for i=1:numel(featdef)
    w = t*featdef{i};
    e{i} = clevalf(t*w,classf,featsizes,[],n,s*w,testfun);
    if numel(classf) > 1
      e{i}.title = [getname(featdef{i}) ' feature curve for ' getname(a)];
    else
      e{i}.names = getname(featdef{i});
    end
  end
  
return
	
