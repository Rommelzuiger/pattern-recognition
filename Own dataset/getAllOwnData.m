function out = getAllOwnData(sequence_range)

    % Define the sequence range
    if (nargin == 0)
        sequence_range = 1:7;
    end

    % First, get all data
    data = {};
    labels = {};
    count_objects = 0;
    for i = sequence_range
        disp(['Starting with loading sequence ' num2str(i)]);
        data_sheets{i} = cropCells(scanDigits(i));
        labels{i} = getLabelsOfSequence(i);
        count_objects = count_objects + length(data_sheets{i});
        disp(['Loading of sequence ' num2str(i) ' is finished']);
    end
    
    % Second, put all data in one object
    data = []; %zeros([count_objects, 16*16]);
    goodlabels = [];
    q = 1;
    for i = sequence_range
        for j = 1:length(data_sheets{i})
            object = im_resize(cell2mat(data_sheets{i}(j)), [16, 16]);
            data = [data; object(:)'];
            goodlabels = [goodlabels; ['digit_' num2str(labels{i}(j))]];
        end
        disp(['Data of sequence ' num2str(i) ' is merged into one big data variable']);
    end
    
    disp(['Final action: putting the data into a prdataset']);
    out = prdataset(data, goodlabels);

end

% Get the right labels of a sequence
function labels = getLabelsOfSequence(sequenceNR)
    load labels.mat
    labels = A(:,:,sequenceNR)';
    labels = labels(:);
end