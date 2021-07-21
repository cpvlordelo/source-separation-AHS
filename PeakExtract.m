function [peakMag, peakFreq, ObSpecMagSmoothed] = PeakExtract(ObSpecMag, FFTfrequencies, para)
% Peak detection directly from Spectrogram
% Input
%   - ObSpecMag             : magnitude spectrogram (in Db) (Each column with the DFT of a frame)
%   - FFTfrequencies        : frequencies where the FFT have been calculated 
%   - para 
%       - peakTh            : the absolute threshold (peaks above the maximum amplitude minus this value are considered)
%       - peakTh_rel        : the relative threshold (peaks above the smoothed envelope plus this value are considered)
%       - peakTh_freq_min   : the initial frequency to look for peaks (peaks above this frequency in bins are considered)
%       - peakTh_freq_max   : the frequency limit (peaks below this frequency in bins are considered)
%       - movL              : the moving average width
%       - typeSmoothing     : type of smoothing function used. Can either be 'Normal Moving Average' or 'Gaussian Moving Average'
%       - sigma             : standard deviation used to create the gaussian filter
%
% Output
%   - peakMag               : cell on which each position has a column vector with the magnitude of the detected peaks for respective frame
%   - peakFreq              : cell on which each position has a column vector with the frequency bins of the detected peaks for respective frame
%   - ObSpecMagSmoothed     : smoothed magnitude spectrogram (in Db)
%
% Original Author: Zhiyao Duan
% Last modified: 9/4/2009
%
% Modified and adapted by: Carlos Lordelo
% Last modified: Jan/2017

peakTh = para.peakTh;                                                       % global threshold
peakTh_rel = para.peakTh_rel;                                               % local threshold
peakTh_freq_min = para.peakTh_freq_min;                                     % initial frequency to look for peaks (frequency bins)
peakTh_freq_max = para.peakTh_freq_max;                                     % frequency limit (frequency bins)
movL = para.movL;                                                           % moving average width used to calculate the smoothed log amplitude spectrum
typeSmoothing = para.typeSmoothing;
sigma = para.sigma;

if sum(ObSpecMag == 0) == length(ObSpecMag)                                 % if all the elements are zero, then return
    return;
end
SpecLength = size(ObSpecMag, 1);
LastFrame = size(ObSpecMag,2);
th = max(0,max(ObSpecMag)-peakTh);                                                 % the global threshold

% calculate smoothed log amplitude spectrum by moving average
ObSpecMagSmoothed = ObSpecMag;

if strcmp(typeSmoothing,'Gaussian Moving Average'),
    for i = 1: LastFrame,
        ObSpecMagSmoothed(:,i) = GaussianMovingAve(ObSpecMag(:,i), movL, sigma);
    end
else
    for i = 1: LastFrame, 
        ObSpecMagSmoothed(:,i) = MovingAve2(ObSpecMag(:,i), movL);
    end
end
    % ObSpecMagSmoothed = smooth(ObSpecMag, movL, 'moving');                    % "smooth" function is in Curve Fitting toolbox. 
%if(IfPlot) 
%    hold on; plot(ObSpecMagSmoothed, 'r');
%    hold on; plot(ObSpecMagSmoothed + peakTh_rel, 'k');
%end

ObSpecMagDiff = ObSpecMag - ObSpecMagSmoothed;                              % calculate the relative log amplitude spectrum

peakMag = {};
peakFreq = {};
for j = 1: LastFrame,
    peakMagFrame = [];
    peakFreqFrame = [];
    %for i = 2: (size(ObSpecMag,1) - 1)
    for i = (peakTh_freq_min): (peakTh_freq_max)
        if (ObSpecMag(i,j) < th(j)) | (ObSpecMagDiff(i,j) < peakTh_rel)                   % a peak should satisfy the global threshold and local threshold
            continue;
        end
        if (ObSpecMag(i,j) < ObSpecMag(i+1, j)) | (ObSpecMag(i, j) < ObSpecMag(i-1, j))       % a peak should be a local maximum
            continue;
        % end if ObSpecMag(i) == max(ObSpecMag(max(i-localRange,1):min(i+localRange,SpecLength)))
        else    
            peakMagFrame = [peakMagFrame; ObSpecMag(i,j)];
            peakFreqFrame = [peakFreqFrame; FFTfrequencies(i)];
        end
    end
    peakMag{j} = peakMagFrame;
    peakFreq{j} = peakFreqFrame;
end
end
% a column vector, peak=1 means that a peak appears at the first positive frequency bin
%if IfPlot
%    hold on; plot(peak, peakAmp, 'ro');
%end