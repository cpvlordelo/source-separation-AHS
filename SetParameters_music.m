%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Set parameters for AHS Separation %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "%--" before a comment of a parameter indicates that this parameter is important
%
% Author: Zhiyao Duan
% Created: 6/20/2012
% Last modified: 5/31/2013

% Modified and adapted by by Carlos Lordelo


AUDIO_FILE = 'Mixture.wav';

%% Parameters for preprocessing input signal
[~, Fs] = audioread(AUDIO_FILE);
numberSources = 2;                              % number of mixed harmonic sources
frameLen_ms = 90;                               % frame length in miliseconds
frameLen =2^nextpow2(frameLen_ms*10^(-3)*Fs);   % frame length in number of bins. 2048 = 93ms for Fs=22050
%zpf = 4;                                       % zero padding factor
win = hamming(frameLen, 'periodic');            % window function
nfft = 2^nextpow2(frameLen);
noverlap = fix(frameLen/2);                     % overlap 50%
%noverlap = fix(3*frameLen/4);                  % overlap 75%
hop = frameLen - noverlap;                      % hop size
%RMS_th = sqrt(0.075);                       % frame whose RMS lower than this threshold will be considered silent

%% Parameters for peak extraction
peakThreshold = 50;                                % global amplitude threshold (dB)
peakThreshold_rel = 8;                             % local amplitude threshold (dB)
peakThreshold_freq_min = 60;                       % minimum frequency to look for a peak (Hz) 60Hz is approx. the first harmonic of 'C2'
peakThreshold_freq_max = 20000;                    % maximum frequency to look for a peak (Hz)
%movL = round(0.01*nfft);                                         % moving average width (bins)
movL = 9;
typeSmoothing = 'Normal Moving Average';           %  type of smoothing function used. Can either be 'Normal Moving Average' or 'Gaussian Moving Average'
%typeSmoothing = 'Gaussian Moving Average';         %  type of smoothing function used. Can either be 'Normal Moving Average' or 'Gaussian Moving Average'
%sigma = (0.1*movL);
sigma = 5;
%% ---------- F0's Estimation ----------------
maxf0Num = numberSources;                      %-- the number of maximum F0s in each frame, highest possible number of sources in a frame
f0min_midi = note2midinum('C2');               %-- lowest possible frequency of F0 (midi number)
f0max_midi = note2midinum('B7');               %-- highest possible frequency of F0 (midi number)
searchRadius_midi = 0.5;                       % radius around each peak to search for f0s (midi)
f0step_midi = 0.1;                             % F0 search step (midi)


%% Parameters for calculating Harmonic Structures
maxHarm = 20;                     %-- The number of harmonics, i.e. harmonic structure feature dimensionality (20 default used)
Thresh_PeakF0Belong = 0.03;       %-- The limit of interval (in frequency ratio, i.e., fpeak/k*f0) to decide if a peak is a harmonic or not (0.03 for half semi-tone)
normEnergy_dB = 100;              %-- The total energy to normalize the harmonic structure of each F0 (dB)

%% Parameters for clustering the Harmonic Structures to find the AHS and HSI
K_min = 2;                  %-- The minimum number of neighbors to be used by the NK clustering algorithm   
K_max = 70;                 %-- The maximum number of neighbors to be used by the NK clustering algorithm 
P_noise_max = 0.2;            %-- Maximum value to use for the parameter to delete noise points from the dataset (Decreasing this value will delete more noise points from the dataset)
P_false = 2;               %-- Parameter to detect and delete false edges between two clusters. (Deacreasing this value causes the algorithm to delete more false edges)
minClusterDensity = 0.3;   %-- Parameter used to identify the significant clusters among all the clusters found in the clusterization step.