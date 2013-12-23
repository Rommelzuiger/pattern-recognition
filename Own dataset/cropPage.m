function I = cropPage(sequence)

im = imread(['Sequence ' num2str(sequence) '.jpeg']);

img = rgb2gray(im);
[~,xi,yi] = roipoly(img);
xMin = min(xi); yMin = min(yi);
xMax = max(xi); yMax = max(yi);

I = imcrop(img, [xMin, yMin, xMax-xMin, yMax-yMin]);


% transform onto square
U = [xi(1:4), yi(1:4)];
goalX = [xi(1) xi(3)];
goalY = [yi(1) yi(3)];
G = [goalX(1) goalX(2) goalX(2) goalX(1); goalY(1) goalY(1) goalY(2) goalY(2)]';
T = maketform('projective', U, G);
I = imtransform(I,T);

imshow(I);

imwrite(I, ['Sequence ' num2str(sequence) ' cropped.jpeg']);