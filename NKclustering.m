function [ clusterPoints_cell, sizeClusters] = NKclustering(HSData, para)
% Cluster the harmonic structures to form the AHS of each instrument.
% Clustering technique follows reference (54) of the article: NK algorithm
% Input
%   - HSData      : Data which the columns are the estimated harmonic
%                   structures of each F0; frame by frame. The number of rows is the dimension of the structures (fixed in 20) 
%   - para
%       - numberSources   : The number of mixed sources, i. e., the number of clusters to be found
%       - K               : The number of nearest neighbors to be used by each point
%       - P_noise         : Parameter for denoising of the inout data
%       - P_false         : Parameter to delete false edges from the mutual neighborhood graph
%
% Output
%   - clusterPoints_cell     : cell on which each position has the points belonging to the same cluster
%   - sizeClusters           : vector with the size of each cluster found
%
% Author: Carlos Lordelo
% Last modified: Jan/2017

%fprintf('Clustering the harmonic structures...');
%% Loading Parameters
numberSources = para.numberSources;
K = para.K;
P_noise = para.P_noise;
P_false = para.P_false;

n = size(HSData, 2); % The number of points to analyze
d = size(HSData, 1); % The dimension of each data point (fixed in 20, in the article)

%% Lambda Knowledge (Eigenvalues of the Local Covariance Matrix) Calculation 

% We first use the Mutual Neighborhhod Value Matrix of all the points in the dataset
% along with other arguments to calculate the eigenvalues of the covariance matrix only.
[eigenvalues, ~] = CalculateLocalFeatures(HSData, CalculateMNV_matrix(HSData), K);
% We will still denoise input data, so we ignore the others outputs.

%% Denoising of Input Data
%threshNoise = mean(mean(eigenvalues, 2), 1) + P_noise*sqrt(mean(var(eigenvalues,1,2),1));
%HSData(:,(mean(eigenvalues,1) > threshNoise)) = [];   % Denoising

threshNoise = nanmean(sum(eigenvalues, 2)./sum(eigenvalues~=0,2), 1) + P_noise*sqrt(mean(var(eigenvalues,1,2),1));
HSData(:,(sum(eigenvalues,1)./sum(eigenvalues~=0,1) > threshNoise)) = [];   % Denoising

n_old = n;
n = size(HSData, 2); % The number of points to analyze is different now

% Constructing new mutual neighborhood and local features
[eigenvalues, neighborIndex] = CalculateLocalFeatures(HSData,CalculateMNV_matrix(HSData), K);

%% Denoising of lambda knowledge

% Calculating inadaptability a(i) of each point in the denoised input space
a = CalculateInadaptability(eigenvalues, neighborIndex, neighborIndex(1,:));

% Checking if neighborhood has false edges using threshold
threshFalse = mean(a) + P_false*sqrt(var(a,1));


indexFalse = find(a > threshFalse);
[~, order] = sort (a(indexFalse), 2, 'descend');
indexFalse = indexFalse(order);
for index = 1:length(indexFalse)
    i = indexFalse(index);
    
    % We have to calculate the new adaptability cause the original value might have changed on previous iterations
    aFalse = CalculateInadaptability(eigenvalues, neighborIndex(:,i), i);       
    
    % If point 'i' still has false edge the while loop is executed
    while aFalse > threshFalse,
        %fprintf('%i, %f , %i \n', index, aFalse, i); % debug
        %fprintf('%i , %i \n', index); % debug
        neighborIndexFalse = neighborIndex(neighborIndex(:,i) > 0,i);       % if there is zero on neighborIndex, it is because the edge has already been cut off
        neighborFalse = HSData(:,neighborIndexFalse);
        neighborFalseCenter = mean(neighborFalse,2);
        
        % Destroying the edge of the most distant point of the neighborhood
        dist = sum((neighborFalse - repmat(neighborFalseCenter,1,size(neighborFalse,2))).^2, 1);
        [~, aux] = max(dist);
        j = neighborIndexFalse(aux);     % point 'j' is the most distance point from center
        
        % Destroying edge between point 'i'and 'j'
        neighborIndex((neighborIndex(:, i ) == j),i) = 0;
        neighborIndex((neighborIndex(:, j ) == i),j) = 0;
        
        %% Recalculating the eigenvalues for the new neighborhood of point 'i' without false edge
        %newNeighborhood = neighborFalse;
        %newNeighborhood(:,aux) = [];
        newNeighborhood = HSData(:,neighborIndex(neighborIndex(:,i) > 0,i));
        %if isequal(newNeighborhood, []),
        if numel(newNeighborhood) == 0,
           % disp ('Neighborhood has zero points, please verify code or increase P.false');
            newNeighborhood = HSData(:,i);
        end
        newNeighborhoodCenter = mean(newNeighborhood,2);
        L = size(newNeighborhood, 2);  % the size of neighborhood can be different than K+1 if some edges had already been cut off
        covarianceMatrix = (1/L).*((newNeighborhood - repmat(newNeighborhoodCenter,1,L))*(newNeighborhood - repmat(newNeighborhoodCenter, 1, L))');
        [eigenvalues(:,i), ~] = sort(eig(covarianceMatrix), 'descend');
        
        %% Recalculating the eigenvalues for the new neighborhood of point 'j' without false edge
        newNeighborhood = HSData(:,neighborIndex(neighborIndex(:,j) > 0,j));
        %if isequal(newNeighborhood, []),
        if numel(newNeighborhood) == 0,
            newNeighborhood = HSData(:,j);
            %disp ('Neighborhood has zero points, please verify code or increase P.false');
        end
        newNeighborhoodCenter = mean(newNeighborhood,2);
        L = size(newNeighborhood, 2);  % the size of neighborhood can be different than K+1 if some edges had already been cut off
        covarianceMatrix = (1/L).*((newNeighborhood - repmat(newNeighborhoodCenter,1,L))*(newNeighborhood - repmat(newNeighborhoodCenter, 1, L))');
        [eigenvalues(:,j), ~] = sort(eig(covarianceMatrix), 'descend');
        %% Recalculating inadaptability
        aFalse = CalculateInadaptability(eigenvalues, neighborIndex(:,i), i);
        %a(i) = a2(i);
        %a(j) = CalculateInadaptability(eigenvalues, neighborIndex(:,j));
        %a = CalculateInadaptability(eigenvalues, neighborIndex);
        %lambdaHat =  mean(eigenvalues(:,neighborhood),2);
        %a(i) = a2(i);
        %a(j) = a2(j);
        %aFalse = a(i);
%         if (a(i) > thresholdFalse),
%             aFalse = a(i);
%         else
%             [aFalse, i] = max(a);
%         end
%         %    control = control+1;
    end
end
%% Clustering in lambda knowledge embedded space
[cluster, sizeClusters, numberClusters] = ClusterData(neighborIndex);

%% Sorting the clusters by density. From higher to lower.
[sizeClusters, order] = sort(sizeClusters, 'descend');
cluster = cluster(order);

clusterPoints_cell = cell(1,numberClusters);
for i = 1:numberClusters,
      clusterPoints_cell{i} = HSData(:,cluster{i});
end

end        


