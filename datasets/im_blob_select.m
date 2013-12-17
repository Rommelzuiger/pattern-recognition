function out = im_blob_select(in)

I = bwlabel(in);
n = max(I(:));
sz = zeros(n,1);
for i=1:n
   sz = sum(I(:)==i);
end
[dummt,mx] = max(sz);
out = in(I==mx);

