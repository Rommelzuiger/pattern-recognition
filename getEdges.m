function [image, bounds] = getEdges(data)

    if (size(data, 1) == 1)
        image = computeEdges(data);
        [y, x] = find(image == 1);
        bounds = [x y];
    else
        for i = 1:size(data, 1)
            image{i} = computeEdges(data(i));
            [y, x] = find(image{i} == 1);
            bounds{i} = [x y];
        end
    end

end

function bounds = computeEdges(image)
    % Get logical matrix
    I = data2im(image); 
    % data2im is the secret to getting a double NxM matrix from a prdataset

    % Zero pad around edges
    % This is needed because otherwise the contour does not consider pixels
    % touching the edge
    S = size(I);
    I = imclose(I, [[1, 1, 1];[1, 1, 1];[1,1,1]]);
    I2 = [zeros(1, S(2) + 2); zeros(S(1),1) I zeros(S(1),1); zeros(1, S(2) + 2)];

    % Return edges
    % Use canny_old because the new methods don't always return a closed
    % contour
    % bw = edge(I2, 'zerocross');
    dilate = imerode(I2, [[1, 1, 1];[1, 1, 1];[1, 1, 1]]);
    bounds = I2 - dilate;
end
