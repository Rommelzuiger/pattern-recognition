function features = normalizeFeatures(features)
    for i = 1:size(features, 2)
        features(:,i) = features(:,i) / max(abs(features(:,i)));
    end
    
end