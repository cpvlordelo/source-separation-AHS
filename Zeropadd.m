function [ column_padded ] = Zeropadd( column, L )
% Transform column vector into a column vector of a greater size L 
% zero padding the additional terms

column_padded = zeros(L, 1);
column_padded(1:length(column)) = column(1:length(column)); 
end

