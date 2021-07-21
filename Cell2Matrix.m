function [ x_array ] = Cell2Matrix( x_cell )
% This function gets a cell (1 x N) whose columns have different lengths
% and transforms it in an array (M x N). 
%
% M is the maximum length of the columns of x_cell
%
% The indexes that didnt exist in x_cell become zero in x_array

N = length(x_cell);

M = 0;
for i = 1:N
    M = max(M,length(x_cell{i}));
end
M_cell = cellfun(@(x) M, cell(size(x_cell)), 'UniformOutput', false);    % we need to copy M in a cell to use cellfun next
aux_cell = cellfun(@Zeropadd, x_cell, M_cell, 'UniformOutput', false);   
x_array = cell2mat(aux_cell);

end

