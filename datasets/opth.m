function opth(x,H)

%----------------------------------------------
% opth(x,H)
%
% Repeats a kernel smoothing train-test 
% procedure 5 times on a given regression 
% dataset for the range of kernel values 
% an displays an errobar of the MSE vs. H.
%
% INPUT:
%     x: regression dataset
%     H: vector of kernel widths to evaluate
%----------------------------------------------

b = 5;           % number of batches
n = size(x,1);   % number of objects
I = (1:n) - 1;
I = mod(I,b) + 1;
for j=1:length(H),    
	for k = 1:b
		% extract training set:
		xtrain = x(find(I~=k),:);
		% extract test set:
		xtest = x(find(I==k),:);
		% train smoother:
		w = ksmoothr(xtrain,H(j));
		% evaluate kernel smoother
		MSE(k,j) = xtest*w*testr;
	end
end
figure
errorbar(H,mean(MSE,1),std(MSE,0,1));
xlabel('H');
ylabel('MSE');
