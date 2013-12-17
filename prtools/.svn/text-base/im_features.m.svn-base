%IM_FEATURES Image features computed from Matlab image toolbox
%
%		F = IM_FEATURES(A,GRAY,FEATURES)
%
% INPUT
%   A        Dataset with binary object images dataset
%   GRAY     Gray-valued images (matched with A, optional)
%   FEATURES Features to be computed
%
% OUTPUT
%   F        Dataset with computed features
%
% This function is a *replacement* for the function IM_MEASURE that is
% based on the DipLib measure.m. This implementation is using the Matlab
% function REGIONPROPS from the image toolbox. Use HELP REGIONPROPS to
% find out which features are exactly supported.
% In each image of the measurement set GRAY the features given in FEATURES 
% are measured. In A a segmented version of GRAY has to be supplied.
% When no GRAY is supplied, the binary images in A are used.
%
% Use FEATURES = 'all' for computing all features.
%
% SEE ALSO
% DATASETS, DATAFILES, REGIONPROPS

% Copyright: D.M.J.Tax, D.M.J.Tax@prtools.org.
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_features(a,gray,features)

	 
if nargin < 3  features = []; end
if nargin < 2  gray = []; end

if nargin < 1 | isempty(a)
   b = prmapping(mfilename,'fixed',{gray,features});
   b = setname(b,'Image features');
elseif isdataset(a)
   if (nargin < 3 | isempty(features)) & ~isdataset(gray)
      features = gray;
      gray = a;
   end
   if ~isdataset(gray)
      error('Binary and gray images should be both datasets')
   end
   fsize = getfeatsize(a);
   if any(getfeatsize(gray) ~= fsize)
      error('Image structures of binary and gray images should be identical')
   end
   if length(fsize) == 2, fsize = [fsize 1]; end
   if size(a,1) ~= size(gray,1)
      error('Same number of binary and gray images expected')
   end
   if isempty(features)
      features = 'all';
   end
   im = data2im(a(1,:));
   [props,featlab] = getlabfromregionprops(im,im,features);
   out = [];
   binim = data2im(a);
   grim = data2im(gray);
   nim = size(a,1)*fsize(3);
   for i=1:size(a,1)
      out(i,:) = regionprops(a,props);
   end
   b = setdat(a,out);
   b = setfeatsize(b,[length(f),fsize(3)]);
   b = setfeatlab(b,featlab);

elseif isdatafile(a)
   if nargin < 3 | isempty(features) & ~isdatafile(gray)
      features = gray;
      gray = a;
   end
   if ~isdatafile(gray)
      error('Binary and gray images should be both datafiles')
   end
   if isempty(features)
      features = 'all';
   end
   im = data2im(a(1,:));
   [props,featlab] = getlabfromregionprops(im,im,features);
   % ok, do it by hand:
   dat = [];
   next=1; oldidx = next;
   [b,next] = readdatafile(a,next);
   while next>0
      dat(oldidx,:) = getfeatfromregionprops(data2im(b),data2im(b),props);
      oldidx = next;
      [b,next] = readdatafile(a,next);
   end
   dat(oldidx,:) = getfeatfromregionprops(data2im(b),data2im(b),props);
   % create a dataset:
   b = prdataset(dat,getlabels(a));
   b = setfeatlab(b,featlab);
   b = setname(b,getname(a));

elseif isa(a,'double') 
   if (nargin < 3 | isempty(features)) & ~isdataset(gray)
      features = gray;
      gray = a;
   end
   [props,featlab] = getlabfromregionprops(a,gray,features);
   b = getfeatfromregionprops(a,gray,props);
else
   error('Wrong input')
end
	
return

% [props,featlab] = getlabfromregionprops(im,gray,features)
function [props,featlab] = getlabfromregionprops(im,gray,features)

% compute the properties:
r = regionprops(im,gray,features);
props = fieldnames(r);
% we do it complicated like this, because bwprops could have been 'all',
% and by this we now get all possible (black and white) properties

% remove the things that we are not interested in:
rmfeatures = {'PixelList', 'SubarrayIdx', 'ConvexHull', ...
   'ConvexImage', 'Image', 'FilledImage', 'PixelIdxList', 'Extrema'};
props = setdiff(props,rmfeatures);

% generate the feature labels
featlab = {};
nr = 0;
for i=1:length(props)
   nrf = size(getfield(r,props{i}),2);
   if nrf>1 % there are more measurements per single property
      for j=1:nrf
         nr = nr+1;
         featlab{nr} = [props{i},num2str(j)];
      end
   else % one measurement per property
      nr = nr+1;
      featlab{nr} = props{i};
   end
end
featlab = strvcat(featlab);

return

%   f = getfeatfromregionprops(im,gray,props)
function f = getfeatfromregionprops(im,gray,props)

r = regionprops(im,gray,props);
f = [];
for i=1:length(props)
   f = [f getfield(r,props{i})];
end

