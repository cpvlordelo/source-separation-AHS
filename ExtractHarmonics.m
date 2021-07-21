function [specHarmdB, specMagOutdB, peakAmpdB, peakFreq ] = ExtractHarmonics(specMagdB, peakAmpdB, peakFreq, AHS, HSI, estF0s, parameters)
% Re-estimation of the F0s frame by frame, using the AHS and HSI
%
% Input:
%   - specMagdB             : the log-magnitude spectrogram of the mixture
%   - peakAmpdB          	: a matrix where each column has the detected peak's amplitudes of a frame (in dB)
%   - peakFreq          	: a matrix where each column has the detected peak's frequencies of a frame (in Hz)
%   - AHS          	        : the Average Harmonic Structure of one of the sources
%   - HSI                   : the respective instability for the AHS (a scalar)
%   - estF0s                : a vector with the estimated F0 for each frame
%   - parameters            
%         - DFTfrequencies       : the frequencies (in Hz) of each bin in FFT
%         - Thresh_PeakF0Belong  : the limit of interval (in frequency ratio, i.e., fpeak/k*f0) to decide if a peak is a harmonic or not (0.03 for half semi-tone)
%       
% Output:
%   - specHarmdB            : the magnitude spectrum (in dB) of the source related to the respective AHS
%   - specMagOutdB          : the log-magnitude spectrogram of the mixture after the source's harmonics extraction
%   - peakAmpdB          	: a matrix where each column has the detected peak's amplitudes of a frame (in dB) after the extraction
%   - peakFreq          	: a matrix where each column has the detected peak's frequencies of a frame (in Hz) after the extraction
%
% Author: Carlos Lordelo
% Last Modified: Jul/2017   

% Loading parameters
Thresh_PeakF0Belong = parameters.Thresh_PeakF0Belong;
DFTfrequencies = parameters.DFTfrequencies;

%maxPeak = size(peakFreq,1);
estF0s = estF0s(:);

[maxHarm , numberSources] = size(AHS);
specHarmdB = -Inf*ones(size(specMagdB));
specMagOutdB = specMagdB;
%% Checking Error
if numberSources > 1,
    error(' The number of sources to be extracted from signal should be one, please check the size of the AHS and the trueF0s vector');
    exit;
end

%% Calculating AHS total energy

totalEnergyAHS = sum(10.^((AHS(AHS>0))/10));                 

lastFrame = length(estF0s);

%% Starting the loop
for frameNum = 1:lastFrame
    f0 = estF0s(frameNum);
    
    if f0 == 0                                        
        continue;
    end        
    %% Eliminating the zeros from the columns
    nonzeroIndex = find(peakFreq(:,frameNum) > 0);
    peakFreq_frame = peakFreq(nonzeroIndex,frameNum);
    peakAmpdB_frame = peakAmpdB(nonzeroIndex,frameNum);
    
    % If there is less than 5 harmonics on frame it should be ignored 
    if numel(peakFreq_frame) < 5,
        continue;
    end
    
    %% Normalizing the peaks amplitude to have same energy as AHS
    totalEnergyFrame = sum(10.^(peakAmpdB_frame/10));
    peakAmp_norm = (10.^(peakAmpdB_frame/10))*totalEnergyAHS/totalEnergyFrame;
    peakAmp_norm = 10*log10(peakAmp_norm);
    %% Normalizing the AHS amplitude to have same energy as the frame's
    AHS_norm = AHS;
    AHS_norm(AHS_norm>0) = (10.^(AHS_norm(AHS_norm>0)/10))*totalEnergyFrame/totalEnergyAHS;
    AHS_norm(AHS_norm>0) = 10*log10(AHS_norm(AHS_norm>0));
    %% Calculating the harmonic number of a peak given the F0
    Ratio = peakFreq_frame / f0; 
    RRatio = round(Ratio);      
    Dev = Ratio./RRatio;
        
        % Use the relative frequency deviation to decide which peak can belong to which F0
        RelDev = Dev - 1;                % the relative frequency deviation of all peaks from harmonics of F0. 1 = F0/F0
        
        % Calculate the harmonics in the frame
        %harmidx_vector = zeros(1, maxHarm);
        for k = 1:maxHarm
            harmidx = find((abs(RelDev)<Thresh_PeakF0Belong & (RRatio == k)) == 1);   % find all peaks that could belong to this harmonic
            if isempty(harmidx)
                continue;                      % if no peak is associated, ignore this harmonic
            elseif AHS(k) == 0;
                continue;
            else
                [~, idx] = min(abs(RelDev(harmidx)));
                harmidx = harmidx(idx);           % if there are multiple peaks, use the closest one
                
                harmAmp = peakAmpdB_frame(harmidx);
                %% If the harmonic amplitude is greater than one HSI in relation to the AHS we should use the value from the AHS as the harm amplitude
                if abs(AHS(k) - peakAmp_norm(harmidx)) > HSI,
                    harmAmp = AHS_norm(k);
                end
                %% Getting the frequency bin of the harmonics
                [~, fBin] = min(abs(DFTfrequencies - peakFreq_frame(harmidx)));
                specMagOutdB(fBin-1:fBin+1,frameNum) = specMagdB(fBin-1:fBin+1,frameNum) - harmAmp;
                specHarmdB(fBin,frameNum) = harmAmp;  
                
                peakAmpdB_frame(harmidx) = max((peakAmpdB_frame(harmidx) - harmAmp),0);
                if peakAmpdB_frame(harmidx) == 0,
                    peakFreq_frame(harmidx) = 0;
                end
            end
        end
            peakFreq(nonzeroIndex,frameNum) = peakFreq_frame;
            peakAmpdB(nonzeroIndex, frameNum) = peakAmpdB_frame;
            
            %% Normalizing the energy of the extracted harmonic source to be the same as the mixture spectrum
%             totalEnergySpecHarmFrame = sum(10.^(specMagOut(:,frameNum)/10));
%             totalEnergySpecMagMix    = sum(10.^(specMagdB (:,frameNum)/10));
%             specMagOut(:,frameNum) = (10.^(specMagOut(:,frameNum)/10))*totalEnergySpecMagMix/totalEnergySpecHarmFrame;
%             specMagOut(:,frameNum) = 10*log10(specMagOut(:,frameNum));
end
    %specMagdB = specMagdB - specMagOut;
end

