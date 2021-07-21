function [F0Feature, F0Location, F0Track] = CalculateHarmonicStructure(PeakData, PeakAmpData, EstF0, para)
% Calculate F0 harmonic structure feature vectors
% Input
%   - PeakData      : Frequency (Hz) of extracted peaks of the audio file, each column is a frame
%   - PeakAmpData   : Amplitudes (dB) of extracted peaks of the audio file, each column is a frame
%   - EstF0         : Multiple F0 estimation results (Hz), each column is a frame
%   - para
%       - maxHarm               : The number of harmonics, i.e. harmonic structure feature dimensionality (20 default used)
%       - Thresh_PeakF0Belong   : The limit of interval (in frequency ratio, i.e., fpeak/k*f0) to decide if a peak is a harmonic or not (0.03 for half semi-tone)
%       - normEnergy_dB         : The total energy to normalize the harmonic structure of each F0 
% Output
%   - F0Feature     : Harmonic structure of this F0. Each row is a feature vector. Total energy of harmonics of each F0 is normalized to "totalF0Energy_dB" dB.
%   - F0Locations
%       - Col 1     : Frame number of analysis
%       - Col 2     : F0 Frequency (Hz)
%       - Col 3     : Energy of its harmonics (10*log10 (sum of squared linear amplitudes))
%   - F0Track       : F0 track label
%
% Author: Zhiyao Duan
% Created: 4/28/2009
% Last modified: 7/2/2012

% Modified and adapted by Carlos Lordelo
% Last modified: Jan/2017

% fprintf('Calculating harmonic structure vectors...');

% parameters
Thresh_PeakF0Belong = para.Thresh_PeakF0Belong;
maxHarm = para.maxHarm;
normEnergy_dB = para.normEnergy_dB;

[MaxPoly, FrameNum] = size(EstF0);

% Initialize several variables
F0Features = zeros(MaxPoly*FrameNum, maxHarm);                      % feature vectors
F0Locations = zeros(MaxPoly*FrameNum, 3);
F0Tracks = zeros(MaxPoly*FrameNum, 1);                              % Initial F0 tracks

% Start processing each frame
feaNum = 0;                                                         % count the number of feature vectors
for i = 1:FrameNum
    if sum(PeakData(:,i)~=0) == 0                                   % this frame is silent or there is no peak in the spectrum
        continue;
    end
    
    peak = PeakData(:, i);
    peakAmp = PeakAmpData(:, i);
    
                                                    
    index = find((peak>0));
    peak = peak(index);                                            % because the peak number in each frame might be different, so there are some 0s
    peakAmp = peakAmp(index);    
       
    %mpeak = peak;                                                   % peak frequency in MIDI number
    %fpeak = midi2hz(mpeak);                                         % peak frequency in Hz
    
    fpeak = peak;
    
    for j = 1:MaxPoly
        if EstF0(j, i) == 0                                        
            continue;
        end        
        
        %mF0 = EstF0(j, i);                                          % F0 in MIDI number
        %fF0 = midi2hz(mF0);                                         % F0 in Hz
        
        fF0 = EstF0(j,i);
        
        Ratio = fpeak / fF0; 
        RRatio = round(Ratio);      % calculate the harmonic number of peak i given F01
        Dev = Ratio./RRatio;
        
        % Use the relative frequency deviation to decide which peak can belong to which F0
        RelDev = Dev - 1;                % the relative frequency deviation of all peaks from harmonics of F01. 1 = F0/F0
        
        % Calculate the linear harmonic structure
        F0HarmStruct = zeros(1, maxHarm);
        for k = 1:maxHarm
            harmidx = find((abs(RelDev)<Thresh_PeakF0Belong & (RRatio == k)) == 1);   % find all peaks that could belong to this harmonic
            if isempty(harmidx)
                F0HarmStruct(k) = 0;                      % if no peak is associated, use 0. This is much better!
            else                
                [~, idx] = min(abs(RelDev(harmidx)));               % if there are multiple peaks, use the closest one
                F0HarmStruct(k) = 10^(peakAmp(harmidx(idx))/20);   	% ATTENTION! FoHarmStruct IS NOT IN dB!! IT IS IN LINEAR SCALE!! 
            end
        end
        
        % normalize harmonic structure to have the total energy 100 dB for each F0's harmonics
        idx = F0HarmStruct ~= 0;
        tmpLinearHarmonic = F0HarmStruct(idx);
        F0Energy = tmpLinearHarmonic.^2;                 
        totalF0Energy = sum(F0Energy);                          % total energy of the harmonics of this F0
        totalF0Energy_dB = 10*log10(totalF0Energy);             % total energy of the harmonics of this F0 (dB)
        normFactorLin = 10^(normEnergy_dB/10);                  % 1e10 for 100 (dB)
        F0Energy = (F0Energy/totalF0Energy)*normFactorLin;      % normalizing the linear energy for 1e10 if 100 dB
        F0HarmStruct(idx) = 10*log10(F0Energy);                 % normalization of the harmonic structures
        
        % compose the feature vector for each F0
        feaNum = feaNum + 1;
        F0Features(feaNum, :) = F0HarmStruct;
        F0Locations(feaNum, 1) = i;
        F0Locations(feaNum, 2) = EstF0(j, i);
        F0Locations(feaNum, 3) = totalF0Energy_dB;
        F0Tracks(feaNum) = j;                                       % An initial track label for each F0
    end
end

% We collected feaNum features in total
F0Feature = F0Features(1:feaNum, :);
F0Location = F0Locations(1:feaNum, :);
F0Track = F0Tracks(1:feaNum);

%fprintf('\n');
end