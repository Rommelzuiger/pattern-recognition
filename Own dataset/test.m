close all
clear all
prwaitbar off
warning off
prwarning off


tic
a = prnist(0:9, 1:500);
[train_set, test_set] = gendat(a, 0.5);
disp(['Getting NIST. done ' num2str(toc) 'ms'] ); tic

mapping = pcam(prdataset(im_resize(train_set, [16 16])), 0.71); 
disp(['PCA mapping. done ' num2str(toc) 'ms'] ); tic
PCA_dataset = prdataset(im_resize(train_set)) * mapping;
disp(['Remapping training data. done ' num2str(toc) 'ms'] ); tic
classifier = parzenc(PCA_dataset);
disp(['Trained parzen. done ' num2str(toc) 'ms'] ); tic

PCA_testset = prdataset(im_resize(test_set, [16 16])) * mapping;
disp(['Remapping testing data. done ' num2str(toc) 'ms'] ); tic
testc(PCA_testset, classifier);
disp(['Testing classifier. done ' num2str(toc) 'ms'] ); tic


% SEQ = 4;
% cells = scanDigits(SEQ);
% disp(['Identifying cells. done ' num2str(toc) 'ms'] ); tic
% seg = cropCells(cells);
% disp(['Cropping cells. done ' num2str(toc) 'ms'] ); tic
% data = convertSeg(seg,SEQ);
% disp(['Converting to PRdata. done ' num2str(toc) 'ms'] ); tic

load data

PCA_own = data*mapping;
disp(['Remapping own test data. done ' num2str(toc) 'ms'] ); tic
testc(PCA_own, classifier);
disp(['Testing own data. done ' num2str(toc) 'ms'] ); tic

