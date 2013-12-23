function out = convertSeg(seg, seqNR)
    load labels.mat
    labels = A(:,:,seqNR)';
    labels = labels(:);
    
    data = zeros([length(seg), 16*16]);
    goodlabels = [];
    for i=1:length(seg)
        s = im_resize(seg{i}, [16, 16]);
        data(i, :) = s(:);
        goodlabels = [goodlabels; ['digit_' num2str(labels(i))]];
    end
    
    out = prdataset(data, goodlabels);
end