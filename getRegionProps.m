function features = getRegionProps(images)
    features_count = 33;
    features = [length(images), features_count];
    for i = 1:length(images)
        features(i, 1:features_count) = getRegionPropsOfOne(images{i}, features_count)';
    end
    
    % Normalize features
    features = normalizeFeatures(features);
end

function features = getRegionPropsOfOne(image, features_count)
    image = data2im(image); 

    props = regionprops(image, 'all');
    
    features = zeros([features_count, 1]);
    features(1) = props.Area;
    features(2) = props.Centroid(1);
    features(3) = props.Centroid(2);
    features(4) = props.MajorAxisLength;
    features(5) = props.MinorAxisLength;
    features(6) = props.Eccentricity;
    features(7) = props.Orientation;
    features(8) = props.ConvexArea;
    features(9) = props.FilledArea;
    features(10) = props.EulerNumber;
    features(11) = props.EquivDiameter;
    features(12) = props.Solidity;
    features(13) = props.Extent;
    features(14) = props.Perimeter;
    features(15) = props.BoundingBox(3);
    features(16) = props.BoundingBox(4);
    features(17:32) = props.Extrema(:);
    features(33) = length(props.ConvexHull);

end