function [ output, normFactor ] = normalizeTotalEnergy(input, normEnergydB)
% Normalize column by column of the input using total energy of each column
% as 'normEnergydB' in dB. The function ignores the values of 0 dB
%
% Input:
%   - input             : each column is the unormalized data. (in dB)
%   - normEnergydB      : total energy of each column after normalization
%   (in dB)
%         
% Output:
%   - output            : normalized data (in db)
%   - normFactor        : linear normalization factor (final energy / initial energy)
%
% Author: Carlos Lordelo
% Last modified: February 2017

normFactor = ones(1,size(input,2));
normFactorLin = 10^(normEnergydB/10); 
energyInput = 10.^(input/10);
for i = 1:size(input,2),
    idx = find(input(:,i) ~= 0);
    totalEnergy = sum(energyInput(idx,i),1);                            
    normFactor(i) = normFactorLin/totalEnergy;
    energyInput(:,i) = energyInput(:,i)*normFactor(i);                    
    input(idx,i) = 10*log10(energyInput(idx,i));                       
end
output = input;
end

