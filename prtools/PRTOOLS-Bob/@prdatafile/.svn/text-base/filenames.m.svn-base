%FILENAMES Get filenames of datafile
%
%		[NAMES,DIRS,ROOT] = FILENAMES(A)
%
% INPUT
%   A         DATAFILE
%
% OUTPUT
%   NAMES     Names of the files in which the objects of A are stored.
%   DIRS      Names of the directories in which the objects of A are
%             stored.   
%   ROOT      Rootdirectory of the datafile
%
% DESCRIPTION
% This routine facilitates the retrieval of the files where the objects of
% A are stored
%
% SEE ALSO
% DATAFILES

function [names,dirs,root] = filenames(a)

	if isempty(a.rootpath)
			root = pwd;
		else
			root = a.rootpath;
	end

	m = size(a,1);
	findex = getident(a,'file_index');
  type = gettype(a);
	names = cell(1,m);
	dirs  = cell(1,m);
	
	for j=1:m
		names{j} = a.files{2}{findex(j,1)}(findex(j,2),:);
		dirs{j}  = a.files{1}{findex(j,1)};
	end
	names= char(names);
	dirs = char(dirs);