%DISEX_PROTSELFD Example of forward prototype selection
%
% This example shows the use of PROTSELFD for a greedy forward
% selection of prototypes from a square dissimilarity matrix in
% order to optimize the representation set.
%
% The final plot shown is a 'feature curve'. This is the error
% as a function of the number of features used (here the number
% of prototypes). The error measure is the mean classification 
% error of the 1-NN rule using the given dissimilarities as 
% distances (KNNDC([],1))

d = readchicken(20,45);  % read dissimilarity dataset
w = protselfd(a);        % forward feature selection
n = size(w,2);           % max number of selected prototypes
                         % random prototype (feature) ranking
v = featsel(size(d,2),randperm(n));

K = [1 2 3 5 7 10 15 20 30 50 70 100 150 200 300 500 700 1000];
K = [K(K<n) n];          % dimensionalities to be checked

% In the next step the feature curve is build. Note that features
% here are prototypes (representation objects). 
% knndc([],1) is used for classification, i.e. use the values in d 
% (or d*w) as distances in the 1-NN rule. 
% Testd is used for evaluation as it is resistant against missing 
% classes (classes not available in the representation set).
% In 10 repititions 50% of the data is used for training and 50%
% for testing. Note that final performances are biased as all data
% is used in this example for prototype selection.

ew = clevalf(d*w,knndc([],1),K,0.5,10,[],testd);
ew.names = 'Forward selection';
ev = clevalf(d*v,knndc([],1),K,0.5,10,[],testd);
ev.names = 'Random selection';
ev.title = getname(d);
ev.xlabel = 'Size of Representation Set';
plote({ew,ev})