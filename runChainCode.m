function runChainCode(dataset, portion, size_images)

    if nargin < 3
        % Default size of images = no resize
        size_images = 100;
        if nargin < 2
            % Default portion = 0.5
            portion = 0.5;
        end
    end
    
    % Define the resize method
    if (size_images <= 0)
        resize_method = 'preserve';
    else
        resize_method = 'nearest';
    end
    
    % Create the features and train- and testset
    dataset = im_resize(dataset, [size_images, size_images], resize_method);
    images = getEdges(dataset);
    features = getChainCodeHist(images);
    dataset = prdataset(features, getlabels(dataset));
    [train_set, test_set] = gendat(dataset, portion);
    
    % Run experiments
    classifiers = {knnc, ldc, parzenc, qdc, nmc, fisherc};
    mappings = train_set*classifiers;
    testc(test_set, mappings)
    

end