function dataset = applyPCA(dataset, portion)
    a = pcam(dataset, portion);
    dataset = dataset * a;
  end