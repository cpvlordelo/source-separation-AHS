function [estimatedF0s, maxLikelihood, BIC] = EstimateF0s(peakFreq_frame, parameters)%, debug)
% Estimate F0s in a frame using Maximum Likelihood + BIC
%
% Input:
%   - peakFreq_frame             : a column vector with the frame's detected peaks frequencies (in Hz)
%   - parameters            
%         - maxf0Num             : the number of maximum F0s in each frame
%         - f0min_midi           : lowest possible frequency of a F0 (midi number)
%         - f0max_midi           : highest possible frequency of a F0 (midi number)
%         - f0step_midi          : F0 search step (midi number) [default is 0,1, i.e., 1/10 of a semi-tone]
%         - searchRadius_midi    : radius around each peak to search for F0s. (midi)
% 
% Output:
%   - Estimatedf0s     : the estimated F0s (Hz) in this frame 
%   - maxLikelihood    : the likelihood for the estimated F0s
%   - BIC              : Bayesian Information Criterion
%
% Author: Carlos Lordelo
% Last modified: Dez/2016

sigma_1 = 0.03;                   % standard deviation for likelihood. Is set to 0.03 to represent half of the semitone range
C_1 = 1/(sqrt(2*pi*sigma_1^2));   % normalization factor for the the gaussian likelihood function (the integral must be equal to unity)

peakFreq_frame = peakFreq_frame(:);       % making sure peakFreq_frame is a column vector
peakFreq_frame = peakFreq_frame(peakFreq_frame > 0);
peakFreq_frame = sort(peakFreq_frame, 'ascend');
if numel(peakFreq_frame) == 0;
    estimatedF0s = 0;
    maxLikelihood = -1;
    BIC = -1;
    return;
end

%% Loading Parameters
maxf0Num = parameters.maxf0Num;
f0min_midi = parameters.f0min_midi;
f0max_midi = parameters.f0max_midi;
f0step_midi = parameters.f0step_midi;
searchRadius_midi = parameters.searchRadius_midi;

%% Initially we should suppose there is only one source in the frame
numberF0s = 1;

%f0initial_midi = max(f0min_midi, hz2midi(peakFreq_frame(1)/2^(1/6)));
%f0vector_midi = [f0initial_midi:f0step_midi:f0max_midi];      % all possibles f0's in midi
%f0vector_hz = midi2hz(f0vector_midi) ;                        % all possibles f0's in Hz
peakFreq_midi = hz2midi(peakFreq_frame);
peakFreq_midi = peakFreq_midi((peakFreq_midi >= f0min_midi) & (peakFreq_midi <= f0max_midi));
if numel(peakFreq_midi) == 0;
    estimatedF0s = 0;
    maxLikelihood = -1;
    BIC = -1;
    return;
end

f0vector_midi = [];
for i = 1:length(peakFreq_midi),
    possiblef0 = peakFreq_midi(i);
    neighborhood = [possiblef0-searchRadius_midi:f0step_midi:possiblef0+searchRadius_midi];
    f0vector_midi = [f0vector_midi,neighborhood];
end
f0vector_hz = midi2hz(sort(f0vector_midi, 'ascend'));

peakFreq_frame_matrix = repmat(peakFreq_frame, 1, length(f0vector_hz));
f0vector_hz_matrix = repmat(f0vector_hz, length(peakFreq_frame), 1);

% Each column of 'freqDeviation' is the equation (10) of the article applied in a possible f0 , i.e. : freqDeviation(:,i) = deviation(peakFreq, f0_vector_hz(i))
freqDeviation = (peakFreq_frame_matrix./f0vector_hz_matrix - round(peakFreq_frame_matrix./f0vector_hz_matrix)) ./ round(peakFreq_frame_matrix./f0vector_hz_matrix);
freqDeviation = freqDeviation.^2;
% See equation (8) of the article
likelihood = C_1*exp(-freqDeviation/(2*sigma_1^2));
likelihood = prod(likelihood,1);    

% Maximum Likelihood Calculation
[maxLikelihood, I] = max(likelihood);

%% The best estimation for a f0 given there is only one source in this frame
estimatedF0s = f0vector_hz(I);

%% Increasing the number of sources in the frame
numberF0s = numberF0s + 1;
while numberF0s <= maxf0Num;
    if numel(peakFreq_frame) < numberF0s,
        numberF0s = numberF0s + 1; 
    else
        freqDeviation = min(freqDeviation, repmat(freqDeviation(:,I), 1, size(freqDeviation ,2) ) );
    
        % This is to exclude the previously founded f0 from the calculations
        %newf0initial_hz = max(f0vector_hz(I), peakFreq_frame(numberF0s)/2^(1/6));
        %indexInitial = find(f0vector_hz >= newf0initial_hz);
        
        %f0vector_hz(1:indexInitial) = []; 
        %freqDeviation(:,1:indexInitial) = [];
        
        % Discover index from peakFreq closest to f0
        [~, index] = min (abs(peakFreq_frame - f0vector_hz(I)));
        
        % This is to exclude the previously founded f0 from the next calculations
        f0vector_hz((index-1)*length(neighborhood)+1:index*length(neighborhood)) = -1;
        freqDeviation((index-1)*length(neighborhood)+1:index*length(neighborhood)) = +inf;
        
        likelihood = C_1*exp(-freqDeviation/(2*sigma_1^2));
        likelihood = prod(likelihood,1);
        [likelihood, I] = max(likelihood);
    
        maxLikelihood = [maxLikelihood; likelihood];
        estimatedF0s = [estimatedF0s; f0vector_hz(I)];
    
        numberF0s = numberF0s + 1;
    end
end

%% Calculating BIC to get the proper number of sources in this frame
BIC = log(maxLikelihood) - (1/2)* [1:maxf0Num]'.* log(length(peakFreq_frame));

%% Getting the correct number of sources by maximizing BIC
[~, numberSources ] = max(BIC);

%% Returning the correct number of estimated f0's in this frame
estimatedF0s = estimatedF0s(1:numberSources);
maxLikelihood = maxLikelihood(1:numberSources);

% Excluding F0 = -1
maxLikelihood = maxLikelihood(estimatedF0s ~= -1);
estimatedF0s = estimatedF0s(estimatedF0s ~= -1);
end

