function [AHS, HSI_vector, HSI] = CalculateAHS(HSdataClusters_cell,numberSources)
% This function generates the Average Harmonic Structure (AHS) and the
% Harmonic Structure Instability (HSI) of each source given the previously
% clustered harmonic structures
% 
%
% Input:
%   - HSdataClusters_cell      : a cell where each position has the Harmonics Structures correctly linked to one of the sources 
%   - numberSources            : the ideal number of sources (instruments)
% 
% Output:
%   - AHS                      : a matrix where each column has the AHS of one of the sources
%   - HSI_vector               : the instability of each harmonic in the AHS
%   - HSI                      : the final instability of the AHS, considering only the harmonics different than zero
%
% Author: Carlos Lordelo
% Last Modified: Jan/2017

numberClusters = length(HSdataClusters_cell);
d = size(HSdataClusters_cell{1},1);

 % Initializing some variables
 AHS = zeros(d,numberSources);
 HSI_vector = zeros(d,numberSources);
 HSI = zeros(1,numberSources);
 for j = 1:numberSources,
      clusterData = HSdataClusters_cell{j};
      %L_max = numel(clusterData(1,clusterData(1,:) > 0));
      L_max = sum(clusterData(1,:) > 0); 
      L_min = 0.3*L_max;
      ahs = zeros(d,1);
      hsi = zeros(d,1);
      for i = 1:d
          %L_i = numel(clusterData(i,clusterData(i,:) > 0));
          L_i = sum(clusterData(i,:) > 0);
          if  (L_i >= L_min),
            ahs(i) = sum(clusterData(i,:))/L_i;
            hsi(i) = sum((clusterData(i,clusterData(i,:) > 0) - ahs(i)).^2)/L_i;
          end
      end
      AHS(:,j) = ahs;
      HSI_vector(:,j) = hsi;
      HSI(j) = sum(hsi)./sum(hsi ~= 0);
 end
HSI = sqrt(HSI);
HSI_vector = sqrt(HSI_vector);

%% Ignoring possibly noisy peaks and renormalizing the AHS Models found 

% After learning the AHS Model of each source we should make sure that the
% harmonics were not caused by noise, deleting the values which are less
% than the 'peakThreshold' and we should re-normalize them

% Initializing some variables
% AHS_norm = AHS;
% HSI_vector_norm = HSI_vector;
% HSI_energy_norm = 10.^(HSI_vector/10);
% HSI_norm = zeros(1,numberSources);
%  for i = 1:size(AHS,2)
%      index = find(AHS(:,i) < max(AHS(:,i)) - peakThreshold); % We use the same threshold we used when calculating the peaks
%      AHS_norm(index,i) = 0;
%      HSI_vector_norm(index,i) = 0;
%      
%      % Re-normalizing the AHS found so it has 'normEnergy_dB' dBs as total Energy
%      [AHS_norm(:,i), normFactor] = normalizeTotalEnergy(AHS_norm(:,i), normEnergy_dB);
%      
%      % Re-normalizing the HSI_vector
%      HSI_energy_norm(: ,i) = HSI_energy_norm(:,i)*normFactor^2;                                            % If the normalized AHS increases by 'normFactor', the normalized 
%      HSI_vector_norm(HSI_vector_norm(:,i) > 0,i) = 10*log10(HSI_energy_norm(HSI_vector_norm(:,i) > 0,i));  % variance will increase by 'normFactor' squared  
%      
%      HSI_norm(i) = mean(HSI_vector_norm(HSI_vector_norm(:,i)>0,i),1);
%      %HSI(i) = mean(HSI_vector(:,i),1);
%  end

end

