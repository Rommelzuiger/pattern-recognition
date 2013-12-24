function runChainCode(original_dataset)
    
    normalized_lengths = [7 9 11 13 15 18 22];
    
    for i = 1:length(normalized_lengths)
        for j = 1:length(directions)
            
            % Display the parameters
            disp(['Normalized length: ' num2str(normalized_lengths(i))]);% ', number of directions: ' num2str(directions(j))]);
            
            % Create the features and train- and testset
            features = getChainCode(original_dataset, normalized_lengths(i), directions(j));
            dataset = prdataset(features, getlabels(original_dataset));
            [train_set, test_set] = gendat(dataset, 0.5);

            % Run experiments
            classifiers = {knnc, ldc, parzenc, qdc, nmc, fisherc};
            mappings = train_set*classifiers;
            testc(test_set, mappings)
            
        end
    end
    

end