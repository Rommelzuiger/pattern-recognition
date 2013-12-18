function out = generateSampleNumbers(N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
out = zeros(6,32, N);
for i = 1:N
    R = randi(10,22)-1;
    R = R(1:6,:);
    out(:,:,i) = [repmat([1 2 3 4 5 6 7 8 9 0], 6,1) R];
end
end

