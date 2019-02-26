% Afroditi Eleftheriou
% 15.02.2018
% A simple function to down sample data, averaging every n values
% downsampling(x,n)
% n = number of values to be averaged
% x = variable to be undersampled

function[x_new]=downsampling(x,n)
mean_n = 0; 
j = 0; %index of downsampled variable 
for i = 1:length(x) % x is the variable I want to downsample
    mean_n = mean_n + x(i);
    if mod(i,n) == 0 
        j = j + 1;
        x_new(j) = mean_n / n; % downsampled variable
        mean_n = 0;
    end
end
end