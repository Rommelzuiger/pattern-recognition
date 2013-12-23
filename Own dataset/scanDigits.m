function out = scanDigits(sequence)

if exist(['Sequence ' num2str(sequence) ' cropped.jpeg'], 'file')
    I = imread(['Sequence ' num2str(sequence) ' cropped.jpeg']);
else
    I = cropPage(sequence);
end

deltaW = ceil(size(I,2) / 16);
colTick = (0:15)*deltaW + 1;

deltaH = ceil(size(I,1) / 12);
rowTick = (0:11)*deltaH + 1;


I = [I 128 * ones(size(I,1), 2*deltaW)];
I = [I; 128 * ones(2*deltaH, size(I,2))];
    
cells = cell(12*16,1);
k = 1;
D = I;

tic
for i = 1:12
    y = rowTick(i);
    for j = 1:16
        x = colTick(j);
        mask = zeros(size(I));
        mask(y:(y+deltaH), x:(x+deltaW)) = 1;        
        if mod(j,2) == mod(i,2)
            D(mask==1) = 256 - I(mask==1);
        end

        I2 = imcrop(I, [x, y, deltaW, deltaH]);
        cells(k) = mat2cell(I2);
        k = k+1;
    end    
end
toc
imshow(D);

out = cells;


