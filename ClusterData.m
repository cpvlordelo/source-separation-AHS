function [ clusters, sizeClusters, numberClusters] = ClusterData(neighborIndex)
% Linking neighbors to form clusters 
%
% Input
%   - neighborIndex       : a matrix with the index of the points that are neighbors of each point in the dataset. The column 
%                           'i' has the indexes of the points who are the neighbors of point 'i'. Remember that each column 'i'
%                           starts with 'i' itself to represent the first neighbor in the neighborhood. If 'neighborIndex(i,j) = 0',
%                           it means that the edge between the neighbors 'j' and 'i' has been cut off.
%

% Output.
%   - clusters            : a cell where each position has the index of points that form a cluster in the dataset
%   - sizeClusters        : the size of each cluster found by the algorithm (a row vector) 
%   - numberClusters      : the number of clusters found by the algorithm (a scalar)
%
% Author: Carlos Lordelo
% Last modified: Fev/2017

n = size(neighborIndex,2);     % Number of points
neighborIndex(1,:) = [1:n];    % Making sure the first neighbor of each point is the proper point itself
K = size(neighborIndex,1) - 1; % Number of neighbors of each point. Subtract 1 is necessary cause first neighbor is the point itself 

%% First we create 'n' clusters with size K+1; One cluster for each point with its neighbors 
clusters = num2cell(neighborIndex,1);
numberClusters = n;
i = 1;

% This loop gets the clusters that have some neighbors in common and link them together
while (i <= numberClusters),
    neighborhood = clusters{i};
    neighborhood = neighborhood(neighborhood > 0);
    sizeNeighborhood = length(neighborhood);
    index = 1;
    while(index <= sizeNeighborhood),
        cluster2link = clusters{neighborhood(index)};
        cluster2link = cluster2link(cluster2link > 0);
        neighborhood = [neighborhood; cluster2link];     % the cluster with this neighbor is linked to the actual cluster
        clusters{neighborhood(index)} = zeros(K+1, 1);   % the cluster that was linked is now replenished with zeros
        index = index + 1;
        sizeNeighborhood = length(neighborhood);
    end
    clusters{i} = unique(neighborhood);
    i = i+1;
end
%% At this point we have linked some clusters together, but the cell 'clusters' still has size 'n'
%  but with a vector full of zeros on the clusters that have already been linked 

% This loop deletes from the cell 'clusters' the clusters that are all zeros
 clusterReal = {};
 j = 1;
 for i = 1:numberClusters,
     aux = clusters{i};
     aux = aux(aux>0);
     if ~isequal(aux, zeros(0,1)),
         clusterReal{j} = aux;
         j = j+1;
     end 
 end
 numberClusters = j - 1;    % the correct number of clusters after the link.  'numberCluster' < 'n' now 
 clusters = clusterReal;
 
 %% Calculating the size of each cluster
 sizeClusters = zeros(1,numberClusters);
 for i = 1:numberClusters,
     sizeClusters(i) = length(clusters{i});
 end

 %% There is still some points that have been put on two or more clusters. So, we should detect these points,
  % attach them only to the cluster with higher density and delete from the others 

  clusters_mat = Cell2Matrix(clusters);     % converting from cell to matrix zeropadding each column when necessary
for i = 1:n,
    index = find(ismember(clusters_mat, i));
    if length(index) > 1,
        [~,col] = ind2sub(size(clusters_mat), index);
        [~, I] = max(sizeClusters(col));
        clusters_mat(index) = 0;
        clusters_mat(index(I)) = i;
    end
end

clusters = num2cell(clusters_mat,1);     % converting back to cell, but with some zeros added

% Deleting the zeros from the clusters
for i = 1: numberClusters,
    clusters{i} = clusters{i}(clusters{i}>0);
end

% clusterReal = {};
% j = 1;
% for i = 1:numberClusters,
%     aux = clusters_mat(:,i);
%     aux = aux(aux>0);
%     if ~isequal(aux, zeros(0,1)),
%         clusterReal{j} = aux;
%         j = j+1;
%     end 
% end
% numberClusters = j - 1;  
% clusters = clusterReal;

%% Recalculating the size of each cluster
sizeClusters = zeros(1,numberClusters);
for i = 1:numberClusters,
    sizeClusters(i) = length(clusters{i});
    if sizeClusters(i) == 0,         % Just to make sure there is not a cluster with size zero
        clusters(i) = [];
        sizeClusters(i) = [];
        numberClusters = numberClusters - 1;
    end
end
end
