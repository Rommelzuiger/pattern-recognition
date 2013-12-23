function out = cropCells(cells)

N = ceil(sqrt(length(cells)));
k = 1;
% figure
out = cells(size(cells),1);
for i = 1:length(cells)
    mat = cells{i};
    % get 10% chunks
    a1 = ceil(size(mat,1)/10);
    a2 = ceil(size(mat,2)/10);
    % remove example number
    mat = imcrop(mat,[0 4*a1 size(mat,2), size(mat,1)-4*a1]);
    % remove bottom whitespace
    mat = imcrop(mat,[0 0, (size(mat,2)-a2), size(mat,2)]);

    matg = 256-mat;
    patch = true(5); 
    
    % maximum filter
    matg = ordfilt2(matg,numel(patch),patch);
    matg = matg > 10;
    
    % remove boundary lines
    % remove the boundary on the left, rotate the image, repeat 4 times.
    for i2=1:4
        col = 1;
        % sometimes there is a line of zero value pixels before the
        % boundary line starts, so we need a look ahead of 5 columns to
        % check if we should stop.
        while   sum(matg(:,col)) > 0.5 * size(matg,1) ...
                || (col < 10 && sum(matg(:,10)) > 0.5 * size(matg,1))
            matg(:,col) = 0;
            col = col + 1;
        end
        matg = imrotate(matg, 90);
    end
    
    % remove areas that are too small
    %L = bwlabel(matg);
    L = bwareaopen(matg, 300);
    L = bwlabel(L);
    
    % remove zero rows and columns
    L( ~any(L,2), : ) = [];  %rows
    L( :, ~any(L,1) ) = [];  %columns
    
%     subplot(N,N, k);
%     imshow(L);
    k = k+1;
    out(i) = mat2cell(L);
end