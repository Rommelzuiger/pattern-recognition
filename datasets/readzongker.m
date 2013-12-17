% if ispc
% 	load s:/pr/home/ela/data/digits_various/dist_digit_zongker.mat
% else
% 	load /data/pr/home/ela/data/digits_various/dist_digit_zongker.mat
% end

function d = readzongker

[path,fil] = fileparts(mfilename);
s = load(fullfile(path,'dist_digit_zongker.mat'));
lab = s.lab;
d = s.d;
d = sqrt(1-d);
d(1:2001:end) = 0;
%d = (d+d')/2;
d = max(d,d');
d = prdataset(d,lab);
d = setfeatlab(d,lab);


