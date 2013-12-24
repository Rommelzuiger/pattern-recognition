function features = normalizeFeatures(features)
    for i = 1:size(features, 2)
        features(:,i) = features(:,i) / max(features(:,i));
    end
    
end