%% Description
% In this very simple code, we load the breast cancer dataset which can be
% obtained at (www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html)
% and run the OBA algorithm on this dataset using the default options. 

%% Demo Script
load breast_cancer.mat %Loads the data X and the labels y
funObj = LogReg(X,y); %Creates the funObj object using the data loaded
lambda = 1/size(X,1); 
x = OBA(LogReg(X,y),lambda);