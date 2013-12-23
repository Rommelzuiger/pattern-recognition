%%   DisTools Table of Contents
%%
%%
%% Characterization of dissimilarity matrices
%% ------------------------------------------
% CHECKEUCL    Check whether a square dissimilarity matrix has a Euclidean behavior
% CHECKTR      Check whether a square dissimilarity matrix obeys triangle inequality
% CORRTR       Correct a square dissimilarity matrix to obey the triangle inequality
% DISCHECK     Dissimilarity matrix check
% DISNORM      Normalization of a dissimilarity matrix
% DISSTAT      Basic statistics of the dissimilarity matrix
% GOFCL        Goodness of clusters/classes separability vs compactness for dissimilarity data
% INTRDIM      Estimate Intrinsic dimension from dissimilarity data
% ISSQUARE     Check whether a matrix is square
% ISSYM        Check whether a matrix is symmetric
% ASYMMETRY    Compute asymmetry of dissimilarity matrix
% NNE          Leave-one-out Nearest Neighbor error on a dissimilarity matrix
% NNERR        Exact expected NN error from a dissimilarity matrix (1)
% NNERROR      Exact expected NN error from a dissimilarity matrix (2)
% VAT          Visual Assessment of cluster Tendency for dissimilarity matrices	
% 
% 
%  Measures
%  -----------------------------------------------
% BINDISTM     Dissimilarity matrix between binary vectors
% BLURDISTM    Blurred Euclidean distance matrix between blobs
%     BLOBBOX  Find box around a binary blob and resample
% CORRDISTM    Distance matrix based on correlations
% COSDISTM     Distance matrix based on inner products
% DPROCRUSTDM  Distance matrix between datasets based on extended Procrustes problem
% EDGEDISTM    Distance matrix between images based on their edges
% EDITDISTM    Edit distance matrix between strings
% EUDISTM      Euclidean distance matrix
% EXPDISTM     Exponential-type of distance matrix
% FLPDISTM     Fast computation of the lp (p > 0)  distance matrix
% HAMDISTM     Hamming distance matrix between binary vectors
% HAUSDM       Hausdorff and modified Hausdorff distance between datasets of image blobs
% JACSIMDISTM  Jaccard-like distance matrix based on similarities
% LPDISTM      l_p (p > 0) distance matrix
% QDISTM       Distance matrix for quantitative variables
% RANKDISTM    Distance matrix between two data sets based on ranking
% SAMDISTM     Distance matrix based on Spectral Angular Mapper (SAM)
% STRKERM      String Kernel Matrix by Lodhi et al
% 
% 
%  Transformations
%  -----------------------------------------------
% DISSIMT      Fixed DISsimilarity-SIMilarity transformation
% KCENTERM     Kernel weighted centering mapping (also for a similarity matrix) 
% MAKESYM      Make a matrix symmetric
% PROXXM       Proximity mapping
% SIGMOID      Element-wise sigmoid tranformation of a matrix
% 
% 
% 
%  Projections
%  -----------------------------------------------
% FASTMAPD     FastMap; inear projection of Euclidean distances
% PE_EM         Pseudo-Euclidean embedding (includes Classical Scaling as a special case)	
% SPHEM        Spherical Embedding
% 
% 
%  Pseudo-Euclidean spaces and indefinite kernels
%  -----------------------------------------------
% IKFD         Indefinite Kernel Fisher discriminant
% IKPCA        Indefinite Kernel PCA
% PE_AFFINE
% SETSIG
% GETSIG
% ISPE_DATASET
% ISPE_EM
% PE_DISTM     Square pseudo-Euclidean distance between two datasets
% PE_KERNELM
% PE_MTIMES
% PE_PARZENC
% PE_KNNC
% PE_NMC
% PE_EM         Pseudo-Euclidean linear embedding
% PSPCA        Pseudo-Euclidean Principal Component Analysis
% 
% 
%  Routines supporting in learning
%  -----------------------------------------------------------------------
% CROSSVALD    Cross-validation error for dissimilarity representations
% CLEVALD
% GENDDAT      Generate random training and test sets for dissimilarity data
% GENREP       Generate a representation set
% GENREPI      Generate indices for representation, learning and testing sets
% SELCDAT      Select Class Subset from a Square Dissimilarity Dataset
% PROTSELFD    Forward prototype selection   
% 
%  Classifiers
%  -----------------------------------------------------------------------
% AUCDLPC      AUC-LP classifier on dissimilarity data
% DLPC         LP-classifier on dissimilarity (proximity) data
% DRSSCC       Dissimilarity-based random subspace combining classifier 
% KNNDC        K-Nearest Neighbor classifier for dissimilarity matrices
% KFD          Kernel Fisher Discriminant
% KSVC         Kernel Support Vector classifier on a kernel matrix
% 	KSVO        Kernel Support Vector Optimizer
% KSVC_NU      Kernel Support Vector classifier on a kernel matrix; nu-version 
% 	KSVO_NU     Kernel Support Vector Optimizer; nu-version 
% MCLASSDC     Multi-Class Dissimilarity-based Classifier from Two-Class Discriminants
% TESTKD       Test k-NN classifier for dissimilarity data
% TQDC         Trade-off Quadratic Discriminant (Regularized Bayes Normal Classifier)
% 
% 
%  Graphs and distances
%  -----------------------------------------------
% DISTGRAPH    Computes distances in a graph
% DMSTSPM      Finds the shortest paths along K minimum spanning trees
% DSPATH       Single shortest path in a (dissimilarity) graph
% DSPATHS      All shortest paths in a (dissimilarity) Graph
% GRAPHPATH    Compute shortest paths in a graph
% KMST         Finds K minimum spanning trees based on a distance matrix
% MSTPLOT      Plot minimum spanning trees
% NHGRAPH      Find a neighborhood graph and its shortest paths
% PLOTGRAPH    Plot a 2D graph
% 
% 
%  Additional
%  -----------------------------------------------
% READCHICKEN  Read chicken data
% PRTBUNKE 	  Load Bunke dissimilarity data to PRTOOLS	

%  Superfluous / outdated but still in CVS
%  -------------------------------------------------
%  KPCA, AUGPSEM, PSEM
%  