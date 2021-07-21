function [eigenvalues, neighborIndex] = CalculateLocalFeatures(data, MNV_matrix, K)
% Calculate local features for all points in the dataset.
%  
% Input
%   - data          : Each column is a point of the dataset 
%   - MNV_matrix    : Mutual Neighborhood Value Matrix giving information of the MNN distances between the points 
%   - K             : Neighborhood size. The number of nearest neighbors to be found for each point [scalar]
% Output
%
%   - eigenvalues          : Each column i has the Eingenvalues of the covariance matrix of the K-neighborhood of point x(i).
%                            Note that eigenvalues(:,i) = [lambda_1; lambda_2; lambda_3; ... lambda_d] such that lambda_1 > 
%                            > lambda_2 > lambda_3 > ... > lambda_d. And d is the dimension of the dataset
%
%   - neighborIndex        : Each column i has the K+1 indexes of variable data to find the neighborhood of x(i) = [x(i) , x(i)_1 , 
%                            , x(i)_2 , ... , x(i)_K], where x(i)_p is the p'th nearest neighbor of x(i). Note that size(neighborIndex) = [(K+1), n]. 
%                            Each column has K+1 as size because it includes the index for finding x(i) itself on the first row of column i
%   
% Author: Carlos Lordelo
% Last modified: Fev/2017

[~, neighborIndex ] = sort(MNV_matrix, 1, 'ascend');
n = size(data,2);      % number of points to analyze
d = size(data,1);      % dimension of thedataset

neighborIndex = neighborIndex((1:K+1),:);

% Preallocating memory
eigenvalues = zeros(d,n);
for i = 1:n,
    
    % Calculating the K nearest neighbors of x_i 
    neighborhood = data(:,neighborIndex(:,i));    
    %neighborhoodCenter(:,i) = mean(neighborhood,2);
    neighborhoodCenter = mean(neighborhood,2);
    
    % Calculating the Covariance Matrix and its eigenvalues for each neighborhood of K elements 
    covarianceMatrix = (1/(K+1)).*((neighborhood - repmat(neighborhoodCenter,1,K+1))*(neighborhood - repmat(neighborhoodCenter,1,K+1))');
    
    % Sorting the eigenvalues such that lambda_1 > lambda_2 > lambda 3 , ... , lambda_d
    [eigenvalues(:,i), ~] = sort(eig(covarianceMatrix), 'descend');
    
end
end

