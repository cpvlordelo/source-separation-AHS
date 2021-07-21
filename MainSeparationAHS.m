
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Main Code for AHS Separation %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Carlos Lordelo
% Last Modified: Jul/2017

close all;
clear all; 

addpath .\MIDIToolbox
%% Getting some important parameters for music analysis
SetParameters_music;

%% Reading Audio File & normalizing mean to zero and RMS value to 1
[wavData, Fs] = ReadWavFile_Mono(AUDIO_FILE);

%% Calculating the spectrogram of the audio signal;
[spec, f, t] = spectrogram(wavData, win, noverlap, nfft, Fs);
%[spec, f, t] = STFT(wavData, frameLen, hop, nfft, Fs);

specMag = abs(spec);
%X_phase = unwrap(angle(X));

% Putting the result in Db and keeping spectrogram as an array
specMagdB = mag2db(specMag);

% Variable representing the spectrogram as a cell, where each position is the spectrum of a frame
specMagdB_cell = mat2cell(specMagdB,size(specMagdB,1), [ones(1,size(specMagdB,2))]);


%% Peak Detection
parameters = {};
parameters.peakTh = peakThreshold;                                                % global threshold for peak detection
parameters.peakTh_rel = peakThreshold_rel;                                        % local threshold for peak detection
parameters.peakTh_freq_min = round(peakThreshold_freq_min*nfft/Fs);               % minimum frequency to search for peaks (frequency bins)
parameters.peakTh_freq_max = min(nfft/2, round(peakThreshold_freq_max*nfft/Fs));  % frequency limit to search for peaks (frequency bins)
parameters.movL = movL;                                                           % width of the moving average filter to smooth spectrum
parameters.typeSmoothing = typeSmoothing;                                         % Type of smoothing function used. Can either be 'Normal Moving Average' or 'Gaussian Moving Average'
parameters.sigma = sigma;                                                         % standard deviation used to create the gaussian filter

% Extracting peaks from SpecMag in dB
[~, peakFreq_cell, specMagdBSmoothed] = PeakExtract(specMagdB,f,parameters);

% Remember SpecMagDbPeak_cell is a cell on which each position has a column
% vector with the magnitude of the detected peaks from its respective frame
% and PeakFreq_cell is a cell on which each position has a column
% vector with the frequencies in Hz of the detected peaks in its respective frame

% Using a parabolic interpolation on the peaks as shown in the article [48]
 f_cell = cellfun(@(x) f, cell(size(peakFreq_cell)), 'UniformOutput', false);    % we need to copy f in a cell to use cellfun next
[peakAmpdB_cell, peakFreq_cell] = cellfun(@PeakInterpolate, specMagdB_cell, peakFreq_cell, f_cell, 'UniformOutput', false);

%% Validation of the previously found peaks

%[peakAmpdB_cell, peakFreq_cell] = cellfun(@PeakValidation, peakAmpdB_cell, peakFreq_cell, specMagdB_cell, f_cell, 'UniformOutput', false);

%% F0 Estimation
parameters = {};
parameters.maxf0Num = maxf0Num;
parameters.f0min_midi = f0min_midi;
parameters.f0max_midi = f0max_midi;
parameters.f0step_midi = f0step_midi;
parameters.searchRadius_midi = searchRadius_midi;

%debug_cell = mat2cell(1:length(peakFreq_cell),[1],ones(1,length(peakFreq_cell)));
parameters_cell = cellfun(@(x) parameters, cell(size(peakFreq_cell)), 'UniformOutput', false);    % we need to copy parameters in a cell to use cellfun next
[estimatedF0s_cell, maxLikelihood_cell, BIC_cell] = cellfun(@EstimateF0s, peakFreq_cell, parameters_cell, 'UniformOutput', false); %debug_cell, 'UniformOutput', false);

% Transforming some cells to matrix variables. 
% Since each cell position has only column vectors with possibly different lengths (diff number of peaks in two frames for example),
% each constructed matrix has "max(length(cell{i})), for i = 1:totalFrameNum" number of rows. The added values in each column are all zeros.

peakAmpdB = Cell2Matrix(peakAmpdB_cell);
peakFreq = Cell2Matrix(peakFreq_cell);
estimatedF0s = Cell2Matrix(estimatedF0s_cell);

%% Features (Harmonic Structure) Extraction
parameters = {};
parameters.maxHarm = maxHarm;
parameters.Thresh_PeakF0Belong = Thresh_PeakF0Belong;
parameters.normEnergy_dB = normEnergy_dB; 

% Calculating the Harmonic Structures of each F0, frame by frame
[harmStruct, F0Location, F0Track] = CalculateHarmonicStructure(peakFreq, peakAmpdB, estimatedF0s, parameters);

% Before learning the AHS Model of each source we should make sure all the harmonic structure have at least 5 proeminent harmonics between the first 20's
totalNumStruct = 0;
for i = 1:size(harmStruct,1)
    s = harmStruct(i,:);
    %s(s < (max(s)-peakThreshold)) = 0;
    if numel((s(s>0)) >= 5) %&& (s(1) > 0),
        totalNumStruct = totalNumStruct + 1;
        HSdata(totalNumStruct,:) = s;
    end
end
HSdata = HSdata';

%   Debug plotting
%    for i = 100:130,
%    figure; plot(f,specMagdB_cell{i}); hold on;
%    stem(peakFreq_cell{i}, peakAmpdB_cell{i}, 'r')
%    plot(f,specMagdBSmoothed(:,i), 'k');
%    end

%% Clustering the Harmonic Structures in order to get an Average Harmonic Structure for each source
parameters = {};
parameters.numberSources = numberSources;
parameters.K = K_min;      
parameters.P_noise = P_noise_max;
parameters.P_false = P_false;

[HSdataClusters_cell, sizeClusters] = NKclustering(HSdata, parameters);

%% Checking if we should change parameters of the clusterization step
parameters.K_min = K_min;
parameters.K_max = K_max;
parameters.minClusterDensity = minClusterDensity;

[HSdataClusters_cell, sizeClusters, K_neighbors, P_noise] = CheckClusterization(HSdata, HSdataClusters_cell, sizeClusters, parameters);

%% Calculating the AHS and HSI for each source and Plotting the results
[AHS, HSI_vector, HSI] = CalculateAHS(HSdataClusters_cell, numberSources);

for i = 1 : numberSources,
figure
stem(AHS(:,i));
end
 
%% Using the AHS of each source to re-estimate the F0s in each frame
parameters = {};
parameters.f0min_midi = f0min_midi;
parameters.f0max_midi = f0max_midi;
parameters.f0step_midi = f0step_midi;
parameters.searchRadius_midi = searchRadius_midi;

[trueF0s, maxLikelihood_matrix] = EstimateTrueF0s2(peakFreq, peakAmpdB, AHS, HSI, estimatedF0s, parameters);

%% Plotting Pianorolls of the sources
parameters = {};
parameters.Fs = Fs;
parameters.frameLen = frameLen;

plotPianoroll(t,trueF0s,parameters)

%% Extracting Harmonics and generating the sources spectrograms
parameters = {};
parameters.DFTfrequencies = f;
parameters.Thresh_PeakF0Belong = Thresh_PeakF0Belong;

specMagHarm_cell = cell(1,numberSources);
source_cell{i} = cell(1,numberSources);

for i = 1:numberSources,
    % Extracting Harmonics of the source
    [specMagHarm_cell{i}, specMagdB, peakAmpdB, peakFreq] = ExtractHarmonics(specMagdB, peakAmpdB, peakFreq, AHS(:,i), HSI(i), trueF0s(i,:), parameters);
    specSource = db2mag(specMagHarm_cell{i}).*exp(j*angle(spec));
    
    % Calculating the Inverse-STFT of the spectrum of the source 
    source_cell{i} = ISTFT(specSource,win,frameLen, hop, 0);
    
    % putting to zero mean
    source_cell{i} = source_cell{i} - mean(source_cell{i});
    % Normalizing the sources 
    source_cell{i} = source_cell{i}/max(abs(source_cell{i}));
end

%% Resynthesis of the rest of the mixture spectrum
specRest = db2mag(specMagdB).*exp(j*angle(spec));
res = ISTFT(specRest, win, frameLen, hop, 0);

% Normalizing the sources and residual signal
res = res/max(abs(res));
wavData = wavData/max(abs(wavData));

%% Putting input and output the same length
if length(wavData) < length(res)
    res = res(1:length(wavData));
    for i = 1:numberSources,
        source_cell{i} = source_cell{i}(1:length(wavData));
    end
else
    res = [res; zeros(length(wavData) - length(res),1)];
    for i = 1:numberSources-1,
        source_cell{i} = [source_cell{i}; zeros(length(wavData) - length(source_cell{i}),1)];
    end
end

%% Saving Results
audiowrite('.\Results\Residual.wav', res,Fs);
nameString = {'1', '2', '3', '4', '5'};
for i = 1:numberSources,
    audiowrite(['.\Results\Source' nameString{i} '.wav'], source_cell{i}, Fs);
end
