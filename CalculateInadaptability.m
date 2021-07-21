function [a] = CalculateInadaptability(eigenvalues, neighborIndex, neighborhoodNumber);
% Calculates the inadaptability of all points in the dataset as equation
% (4) of reference [54] of the article
% Input
%   - eigenvalues          : Each column i has the Eingenvalues of the covariance matrix of the K-neighborhood of point x(i).
%                            Note that eigenvalues(:,i) = [lambda_1; lambda_2; lambda_3; ... lambda_d] such that lambda_1 > 
%                            > lambda_2 > lambda_3 > ... > lambda_d. And d is the dimension of the dataset
%
%   - neighborIndex        : Each column i has the K+1 indexes of variable data to find the neighborhood of x(i) = [x(i) ; x(i)_1 ; 
%                            ; x(i)_2 ; ... ; x(i)_K], where x(i)_p is the p'th nearest neighbor of x(i). Note that size(neighborIndex) = [(K+1), n]. 
%                            Each column has K+1 as size because it includes the index for finding x(i) itself on the first row of column i
%
%   - neighborhoodNumber   : Necessary for calculations of inadaptabilities for only parts of the dataset, i.e., if 'neighborIndex' has different 
%                            length than 'eigenvalues'. This variable has the first neighbor for each column of 'neighborIndex' so it is possible
%                            to control the index of 'eigenvalues' we should calculate 'a' with
% Output
%   - a                    : Inadaptability of the points in the dataset
%   
% Author: Carlos Lordelo
% Last modified: Fev/2017

n = size(neighborIndex,2);
K = size(neighborIndex, 1) - 1;
d = size(eigenvalues, 1);

a = zeros(1,n);
for i = 1:n
    neighborhood = neighborIndex(2:end,i); % we need to begin on 2 because we should take the point x_i itself (first position in each column) out of the calculation of the mean
    neighborhood = neighborhood(neighborhood > 0); % there may be zeros on neighborIndex indicating the edge between point i and j has been cut off
    % neighborhood = neighborhood(neighborhood ~= i);
    % Calculating lambdaHat_ij
    %lambdaHat =  mean(eigenvalues(:,neighborIndex((neighborIndex(2:end,i) > 0),i)),2); %  of the calculation of the mean
    lambda = eigenvalues(:,neighborhood);
    lambdaHat =  sum(lambda,2)./sum(lambda~=0, 2);
    a(i) = nanmean((eigenvalues(:,neighborhoodNumber(i)) ./ lambdaHat), 1);
end

end