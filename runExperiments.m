PCA_dataset = applyPCA(dataset, 0.95);
[train_set, test_set] = gendat(dataset, 0.5);
classifiers = {knnc, ldc, parzenc, qdc, nmc, fisherc};
mappings = train_set*classifiers;
testc(test_set, mappings)