function [MNV] = CalculateMNV_matrix(data)
% Calculate the Mutual Neighborhhod Value Matrix for all points in dataset.
% If x_i is the p'th nearest neighbor of x_j and x_j is the q'th nearest neighbor of x_i,
% MNV_ij = p + q.
%
% Input
%   - data      : Each column is a point of the dataset 
% Output
%   - MNV       : Mutual Neighborhood Value Matrix of the dataset
%   
% Author: Carlos Lordelo
% Last modified: Fev/2017

n = size(data, 2); % The number of points to analyze

MNV= zeros(n,n); % Preallocating memory

for i = 1:n,
    x_i = data(:,i);
    
    % Calculating Euclidean distance of each point
    dist = sqrt(sum((data - repmat(x_i,1,n)).^2, 1));
    
    % Generating the Euclidean Distance Matrix
    %distMatrix = [distMatrix; dist];
    
    % Sort row by distance
    [~, O] = sort(dist, 'ascend');
    
    % Mutual Neighborhood Value Matrix
    MNV(i,O) = [0:(n-1)];
end

MNV = MNV + MNV';

end

