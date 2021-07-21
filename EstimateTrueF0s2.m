function [trueF0s, maxLikelihood_matrix] = EstimateTrueF0s2(peakFreq, peakAmpdB, AHS, HSI, estimatedF0s, parameters)
% Re-estimation of the F0s frame by frame, using the AHS and HSI
%
% Input:
%   - peakFreq          	: a matrix where each column has the detected peak's frequencies of a frame (in Hz)
%   - peakAmpdB          	: a matrix where each column has the detected peak's amplitudes of a frame (in dB)
%   - AHS          	        : a matrix where each column has the AHS of one of the sources
%   - HSI                   : The respective instability for each AHS
%   - estimatedF0s          : a matrix where each row has the first stage estimated F0s for each source
%   - parameters            
%         - f0min_midi           : lowest possible frequency of a F0 (midi number)
%         - f0max_midi           : highest possible frequency of a F0 (midi number)
%         - f0step_midi          : F0 search step (midi number) [default is 0,2, i.e., 1/5 of a semi-tone]
%         - searchRadius_midi    : radius around each peak to search for F0s. (midi)
% 
% Output:
%   - trueF0s                    : A (numberofSources X numberofFrames) Matrix with the re-estimated F0s (in Hz) frame by frame. 
%                                  The first F0 in each column is calculated with respect to the first instrument in the AHS matrix, the second
%                                  F0 is with respect to the second instrument and so on...
%
%   - maxLikelihood_matrix       : maxLikelihood_matrix(i,j) is the maximum likelihood found for the F0 calculated in 'trueF0s' for the instrument i on frame j  
%
% Author: Carlos Lordelo
% Last Modified: Jul/2017

%% Checking error
if size(peakFreq) ~= size(peakAmpdB),
    error('ERROR: The size of the peak-frequency vector is different from the peak-amplitude vector in a frame.')
end

%% Loading Parameters
f0min_midi = parameters.f0min_midi;
f0max_midi = parameters.f0max_midi;
f0step_midi = parameters.f0step_midi;
searchRadius_midi = parameters.searchRadius_midi;

%% Creating some important variables 
numberSources = size(AHS,2);
numberFrames = size(peakFreq,2);
maxHarm = size(AHS,1);

sigma_1 = 0.03;                       % standard deviation for likelihood given the peaks frequency only. Is set to 0.03 to represent half of the semitone range
%C1 = (sqrt(2*pi*sigma_1^2));%*(1/(sqrt(2*pi*sigma_2^2)));   % normalization factor for the the gaussian likelihood function (the integral must be equal to unity)

%% Initializing output
trueF0s = zeros(numberSources, numberFrames);
maxLikelihood_matrix = ones(numberSources, numberFrames);
%sigma_2 =

%% Making a processing loop for each frame
for i = 1:numberFrames
    peakFreq_frame = peakFreq(:,i);
    peakAmpdB_frame = peakAmpdB(:,i);
    index = find(peakFreq_frame > 0);
    peakFreq_frame = peakFreq_frame(index);        % Eliminating the previously added zeros from the frequency vector of peaks
    peakAmpdB_frame = peakAmpdB_frame(index);      % Eliminating the previously added zeros from the amplitude vector of peaks
    
    % Just Checking if there is less than 5 peaks available in this frame (no possible F0, only noise peaks)
    if numel(peakFreq_frame) < 5;
        continue;
    end
    peakFreq_frame_midi = hz2midi(peakFreq_frame); % Putting the values of the frequencies in midi in order to calculate the boundaries of F0s
    peakFreq_frame_midi = peakFreq_frame_midi((peakFreq_frame_midi >= f0min_midi) & (peakFreq_frame_midi <= f0max_midi));
    
    % Just Checking if there is no peak available in this frame (no possible F0)
    if numel(peakFreq_frame_midi) < 5;
        continue;
    end
    
    %% Creating the possibles F0s for this frame
    f0vector_midi = [];
    for j = 1:length(peakFreq_frame_midi),
        possiblef0 = peakFreq_frame_midi(j);
        neighborhood = [possiblef0-searchRadius_midi:f0step_midi:possiblef0+searchRadius_midi];
        f0vector_midi = [f0vector_midi,neighborhood];               % All the possibles frequencies (in midi) for F0s
    end
    f0vector_hz = midi2hz(sort(f0vector_midi, 'ascend'));           % All the possibles frequencies (in Hz) for F0s
    %f0vector_hz = [f0vector_hz];
    
    %% Creating the part of the likelihood that depends of the frequency of the peaks only (Equation (14) of the article) p(f_i | f_0)
    peakFreq_frame_matrix = repmat(peakFreq_frame, 1, length(f0vector_hz));
    f0vector_hz_matrix = repmat(f0vector_hz, length(peakFreq_frame), 1);
    
    % Each column of 'freqDeviation' is the equation (10) of the article applied in a possible f0 , i.e. : freqDeviation(:,i) = deviation(peakFreq, f0_vector_hz(i))
    freqDeviation = (peakFreq_frame_matrix./f0vector_hz_matrix - round(peakFreq_frame_matrix./f0vector_hz_matrix)) ./ round(peakFreq_frame_matrix./f0vector_hz_matrix);
    % See equation (15) of the article
    freqDeviation = min(freqDeviation.^2, 4*sigma_1^2);
    %freqDeviation = freqDeviation.^2;
    likelihood_freq = exp(-freqDeviation/(2*sigma_1^2))./sqrt(2*pi*sigma_1^2);
    allPeaksLikelihood_freq = prod(likelihood_freq,1);

    %% Analyzing now the part of the likelihood that depends of the AHS models found (Equation (14) of the article) p(A_i | f_i,f_0,AHS)
    
    peakAmpdB_frame_matrix = repmat(peakAmpdB_frame, 1, length(f0vector_hz));
    Ratio = (peakFreq_frame_matrix./f0vector_hz_matrix);      
    RRatio = round(Ratio);                          % RRatio(i,j) = harmonic number of peak i given F0 = f0vector_hz(j)
    
    % At this point we have the matrix RRatio which gives us the harmonic number for each peak. However we can't represent a single 
    % harmonic by more than 1 peak. So we should exclude some peaks from the F0 estimation as they can't be a true harmonic 
    dist = (Ratio - RRatio).^2;
     for harmNumber = 1:maxHarm
         for k = 1:size(RRatio,2)
             aux = RRatio(:,k);
             idx = find(aux == harmNumber);
             if numel(idx) >= 2,                     % There exist more than 1 peak for this harmonic (WE SHOULD KEEP THE CLOSEST ONE ONLY)
                 dev = dist(idx,k);
                 [~, closest] = min(dev);
                 idx(closest) = [];                  % Keep the closest (Do not modify its value)
                 aux(idx) = 0;                       % Ignoring the other peaks
                 RRatio(:,k) = aux;
             end
         end
     end
    
    RRatio(isnan(RRatio)) = 0;
    RRatio(RRatio > maxHarm) = 0;                   % If some peaks are responsible for harmonics higher than 'maxHarm', we should ignore them as well
    RRatio = RRatio + 1;                            % The +1 is required because because RRatio can be zero, so we have to add 1 to shift every position by 1
     
    allPeaksLikelihood_amp = zeros(1,length(f0vector_hz));
    for index = 1:numberSources,
        B = [inf; AHS(:,index)];                   % We have to shift B down and put inf at first position to take into account the shifted version of RRatio
        %sigma_2_vector = [0;HSI(:,index)]';
        sigma_2 = HSI(:,index);
        %sigma_2 = sqrt(sum(sigma_2.^2)/sum(sigma_2 ~= 0));
        
        for f0pos = 1:size(RRatio,2)
            harmonics = RRatio(:,f0pos);
            idx = find(harmonics > 1);
            peak = peakAmpdB_frame;
            peakAHS = peak(idx);           % Part of the peaks from frame that are possible harmonics of f0vector_hz(f0pos)
            
            % Putting the total energy of the important peaks the same as the source AHS
            totalEnergyB = sum(10.^(B(2:end)./10));
            energyPeakAHS = 10.^(peakAHS./10);
            totalEnergyPeakAHS = sum(energyPeakAHS);
            peakAHS_norm = 10*log10((energyPeakAHS./totalEnergyPeakAHS)*totalEnergyB);
            peak(idx) = peakAHS_norm;
            
            ampDeviation = (peak - B(harmonics)).^2;
            ampDeviation = min(ampDeviation, 4*sigma_2.^2);
            likelihood_amp = exp(-ampDeviation./(2*sigma_2.^2))./sqrt(2*pi*sigma_2^2);            
            allPeaksLikelihood_amp(f0pos) = prod(likelihood_amp);
        end
        [maximumLikelihood, I] = max(allPeaksLikelihood_freq.*allPeaksLikelihood_amp);
        bestF0 = f0vector_hz(I);
        trueF0s(index,i) = bestF0;
        maxLikelihood_matrix(index,i) = maximumLikelihood;
    end   
    numberSilentSources = sum(estimatedF0s(:,i) == 0);  % getting the number of silent sources in this frame
    if (numberSilentSources > 0),
        [~, silentSources] = sort(maxLikelihood_matrix(:,i), 'ascend');     
        trueF0s(silentSources(1:numberSilentSources),i) = 0;                  % put the lowest likelihood sources to zero
    end
end
    %% Using a cascade of length 3 and length 7 median filters in the F0 lines
    trueF0s = medfilt1(medfilt1(trueF0s, 3, [], 2), 7, [], 2);

end