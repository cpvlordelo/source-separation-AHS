function [HSdataClusters_cell, sizeClusters, K_neighbors, P_noise] = CheckClusterization(HSdata, HSdataClusters_cell, sizeClusters, parameters)
% Changes the clusterization parameters in order to make the clusterization step sucessfull.
% Note that the parameter P_false is kept constant.
%
% Input:
%   - HSdata                     : a matrix where each column has the Harmonic Structure found in each frame
%   - HSdataClusters_cell     	 : a cell where each position has the Harmonics Structures previously linked to one of the sources 
%   - sizeClusters          	 : a vector with the number of points in each cluster
%   - parameters            
%         - numberSources        : the ideal number of sources (instruments)
%         - K_max                : the maximun value of neighbors possible to be used by the NK clusterization algorithm
%         - minClusterDensity    : the minimal density of points in a cluster for it being considered significant 
% 
% Output:
%   - HSdataClusters_cell        : a cell where each position has the Harmonics Structures corected linked to its respective source
%   - sizeClusters               : a vector with the final number of points in each cluster
%   - K_neighbors                : the final value of neighbors used in the corrected NK clustering algorithm step
%   - P_noise                    : the final value of the parameter P_noise used in the corrected NK clustering algorithm step
%
% Author: Carlos Lordelo
% Last Modified: Jul/2017

%% Loading some variables
numberSources = parameters.numberSources;
K_max = parameters.K_max;
minClusterDensity = parameters.minClusterDensity;
%numberPoints = size(HSdata,2);

lookForBestK = true;            %Variable used to stop the following while loop
  while (lookForBestK),
      numberClusters = length(sizeClusters);      
      %% This condition is true if we find less clusters than necessary 
      
      % We should return the value of K to K_min and decrease P_noise.
      
      if numberClusters < numberSources,    
          if parameters.P_noise <= 0,
              disp (['Error in clusterization step. Could only find less than ', num2str(numberSources), ' clusters']);
              disp ('Please check if the number of sources is correct or decrease the minimal value');
              disp ('for the density of a significant cluster');
              return;
          end
          parameters.K = parameters.K_min;
          parameters.P_noise = parameters.P_noise - 0.1;
          clc;
          disp (['Decreasing the value of P_noise to ', num2str(parameters.P_noise), '...']);
          [HSdataClusters_cell, sizeClusters] = NKclustering(HSdata, parameters);
      end
      
      %% This condition is true if we find tha correct number of clusters
      
      % We should make sure the most sparse cluster has at least 'minClusterDensity' 
      % of all the total points, otherwise we should decrease P_noise
      
      if (numberClusters == numberSources),
        if (sizeClusters(numberClusters) >= minClusterDensity*sum(sizeClusters))
            lookForBestK = false;
        else
            if parameters.P_noise <= 0,
                disp (['Error in clusterization step. Could only find less than ', num2str(numberSources), ' clusters']);
                disp ('Please check if the number of sources is correct or decrease the minimal value');
                disp ('for the density of a significant cluster');
                return;
            end
            parameters.K = parameters.K_min;
            parameters.P_noise = parameters.P_noise - 0.1;
            clc;
            disp (['Decreasing the value of P_noise to ', num2str(parameters.P_noise), '...']);
            [HSdataClusters_cell, sizeClusters] = NKclustering(HSdata, parameters);
         end
      end
      
      %% This condition is true if we find more clusters than necessary 
      
      % Remember we might have just found some clusters whose sizes are 
      % insignificant and negligible.
      
      % A cluster can be considered negligible if its density is less than
      % 'minClusterDenity'
 
      % we should also make sure the most sparse cluster has at least 'minClusterDensity'
      % of all the total points, otherwise we should still increase the number of neighbors
      
      if (numberClusters > numberSources)
        if (sizeClusters(numberSources) >= minClusterDensity*sum(sizeClusters)) 
            if (sizeClusters(numberSources + 1) < minClusterDensity*sum(sizeClusters))
                lookForBestK = false;        
            else
                switch (parameters.K == K_max),
                case true
                    if parameters.P_noise <= 0,
                        disp (['Error in clusterization step. Could only find more than ', num2str(numberSources), ' clusters.']);
                        disp ('Please check if the number of sources is correct or increase the minimal value'); 
                        disp ('for the density of a significant cluster');
                        return;
                    end
                    parameters.K = parameters.K_min;
                    parameters.P_noise = parameters.P_noise - 0.1;
                    clc;
                    disp (['Decreasing the value of P_noise to ', num2str(parameters.P_noise), '...']);
                    [HSdataClusters_cell, sizeClusters] = NKclustering(HSdata, parameters);
                case false
                    parameters.K = parameters.K + 1;
                    disp (['Increasing the value of K to ', num2str(parameters.K), '.'])
                    [HSdataClusters_cell, sizeClusters] = NKclustering(HSdata, parameters);
                end
            end
        else
            switch (parameters.K == K_max),
            case true    
                if parameters.P_noise <= 0,
                    disp (['Error in clusterization step. Could only find less than ', num2str(numberSources), ' clusters.']);
                    disp ('Please check if the number of sources is correct or decrease the minimal value'); 
                    disp ('for the density of a significant cluster');
                    return;
                end
                parameters.K = parameters.K_min;
                parameters.P_noise = parameters.P_noise - 0.1;
                clc;
                disp (['Decreasing the value of P_noise to ', num2str(parameters.P_noise), '...']);
                [HSdataClusters_cell, sizeClusters] = NKclustering(HSdata, parameters);
            case false
                parameters.K = parameters.K + 1;
                disp (['Increasing the value of K to ', num2str(parameters.K), '.'])
                [HSdataClusters_cell, sizeClusters] = NKclustering(HSdata, parameters);
            end
         end
      end
  end     
 K_neighbors = parameters.K;
 P_noise = parameters.P_noise;
end

