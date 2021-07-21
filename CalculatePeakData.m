function [PeakData, PeakAmpData, PeakRelAmpData, PeakNum] = CalculatePeakData(specData, FrameIndex, para)
% Calculate Peak Data
% Input:
%   - SpecData          : complex spectrogram
%   - FrameIndex        : indicate of which frame is not a noise frame
%   - para
%       - PeakTh        : the absolute threshold (peaks above the maximum amplitude minus this value are considered)
%       - peakTh_rel    : the relative threshold (peaks above the smoothed envelope plus this value are considered)
%       - peakTh_freq   : the frequency limit (peaks below this frequency in bins are considered)
%       - movL          : the moving average threshold
%       - localRange    : the radius of local range to get the local maximum
% Output:
%   - PeakData          : position data
%   - PeakAmpData       : amplitude data
%   - PeakRelAmpData    : relative amplitude data
% 
% Author: Zhiyao Duan
% Last modified: 5/24/2013
%
% Modified and Adapted by 
% Carlos Lordelo in Jan/2017

fprintf('Detecting peaks in each frame...');

[binL, L] = size(specData);                                                 % [number of frequency bins, frame number]

PeakData = zeros(100, L);                                                   % set the upper bound of the number of peaks to 100
PeakAmpData = zeros(100, L);
PeakRelAmpData = zeros(100, L);
PeakNum = zeros(1, L);

IfPlot = 0;                                                                 % show peak detection results in each frame or not

for FrameNum = 1:L
    if mod(FrameNum, 10) == 0
        fprintf('.');
    end
    if mod(FrameNum, 1000) == 0
        fprintf('\n');
    end
    
    if FrameIndex(FrameNum) == 1
        ObSpecMag = 20*log10(abs(specData(2:binL, FrameNum)));              % log amplitude spectrum
        ObSpecMag(ObSpecMag == -Inf) = min(ObSpecMag(ObSpecMag ~= -Inf));   % set -Inf values to the minimum non -Inf value
        
        [peak, peakNum, SmObSpecMag] = PeakExtract(ObSpecMag, para, IfPlot);% peak extraction
        [peak, peakAmp] = PeakInterpolate(ObSpecMag, peak);                 % peak interpolation
                
        peakRelAmp = peakAmp - SmObSpecMag(round(peak));                    % the relative log amplitude spectrum
        
        PeakData(1:peakNum, FrameNum) = peak';
        PeakAmpData(1:peakNum, FrameNum) = peakAmp';
        PeakRelAmpData(1:peakNum, FrameNum) = peakRelAmp';
        PeakNum(FrameNum) = peakNum;
    end
end
fprintf('\n');