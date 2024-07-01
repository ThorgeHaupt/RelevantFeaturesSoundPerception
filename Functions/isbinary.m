function output = isbinary(x,dim)
% this function returns a bool value if the input vector is binary
output = numel(unique(x)) == 2;
