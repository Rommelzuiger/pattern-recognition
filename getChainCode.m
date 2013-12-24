function features = getChainCode(images, normalized_length)
    
    % Define default normalized_length
    if (nargin < 2)
        normalized_length = 20;
    end
    
    % Get features of each image
    features = [length(images), normalized_length + 8];
    for i = 1:length(images)
        features(i, 1:normalized_length + 8) = getChainCodeOfOne(images{i}, normalized_length)';
    end
    
    % Normalize features
    features = normalizeFeatures(features);
end

function features = getChainCodeOfOne(image, normalized_length)
    
    % Intialization
    features = zeros([normalized_length + 8, 1]);

    % Find trace
    image = data2im(image);
    [x,y] = find(image == 1, 1);
    trace = bwtraceboundary(image, [x y], 'nw', 8, Inf, 'clockwise');
    
    % Find angle between points in the trace
    resized_trace = imresize(trace, [normalized_length 2]);
    for i = 1:normalized_length
       
        % Define the index of the next point 
        j = i + 1;
        if (j == normalized_length + 1)
            j = 1;
        end
        
        % Define the angle of the points
        features(i) = sin((computeAngle(resized_trace(i, 1) - resized_trace(j, 1), resized_trace(i, 2) - resized_trace(j, 2)) / 360) * 2 * pi);
        
    end
    
    % Create histogram
    hist = zeros([8, 1]);
    for i = 1:length(trace)
        
        % Define the index of the next pixel 
        j = i + 1;
        if (j == length(trace) + 1)
            j = 1;
        end
        
        % Top left
        if trace(i,1) == trace(j, 1) - 1 && trace(i, 2) == trace(j, 2) - 1
            hist(1) = hist(1) + 1;
        % Top
        elseif trace(i,1) == trace(j, 1) && trace(i, 2) == trace(j, 2) - 1
            hist(2) = hist(2) + 1;
        % Top right
        elseif trace(i,1) == trace(j, 1) + 1 && trace(i, 2) == trace(j, 2) - 1
            hist(3) = hist(3) + 1;
        % Right
        elseif trace(i,1) == trace(j, 1) + 1 && trace(i, 2) == trace(j, 2)
            hist(4) = hist(4) + 1;
        % Bottom right
        elseif trace(i,1) == trace(j, 1) + 1 && trace(i, 2) == trace(j, 2) + 1
            hist(5) = hist(5) + 1;
        % Bottom
        elseif trace(i,1) == trace(j, 1) && trace(i, 2) == trace(j, 2) + 1
            hist(6) = hist(6) + 1;
        % Bottom left
        elseif trace(i,1) == trace(j, 1) - 1 && trace(i, 2) == trace(j, 2) + 1
            hist(7) = hist(7) + 1;
        % Left
        elseif trace(i,1) == trace(j, 1) - 1 && trace(i, 2) == trace(j, 2)
            hist(8) = hist(8) + 1;
        end
    end
    
    features(normalized_length + 1:normalized_length + 8) = hist(:);
end

% Return the direction of a vector [x y];
function direction = computeDirection(x, y, number_of_directions)
    
    % Define the angle of a part
    angle_of_a_part = 360 / (number_of_directions * 2);
    
    % Compute angle
    angle = computeAngle(x, y);
    
    % Hook to avoid the step between 359.9999 and 0 degrees
    if (angle < angle_of_a_part)
        angle = angle + 360;
    end
    
    % Check the direction
    for direction = 1:number_of_directions
        angle_of_direction = direction * angle_of_a_part * 2;
        if (angle >= angle_of_direction - angle_of_a_part && angle < angle_of_direction + angle_of_a_part)
            return
        end
    end
end


% Return the angle of a vector [x y];
function angle = computeAngle(x, y)
    if atan2(y, x) >= 0
        angle = (180 / pi) * atan2(y, x);
    else
        angle = (180 / pi) * (atan2(y, x) + 2 * pi);
    end
end