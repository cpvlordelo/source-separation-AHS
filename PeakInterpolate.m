function [newSpecMagDbPeak_frame, newPeakFreq_frame] = PeakInterpolate(spectrumMagDb_frame, peakFreq_frame, DFTFrequencies)
% Quadratic interpolation for peak detection
%
% Input
%   - spectrumMagDb_frame      : magnitude spectrum of a frame
%   - peakFreq_frame           : frequency (in Hz) of the frame's previously detected peaks
%   - DFTfrequencies           : vector with frequencies (in Hz) corresponding to the
%                                DFT frequency bins
% Output
%   - newPeakFreq_frame        : frequency (in Hz) of the newly detected peaks of the frame
%   - newSpecMagDbPeak_frame   : amplitudes in dB of the newly interpolated peaks of the frame
% 
% Original Author: Zhiyao Duan
% Created: 2007
%
% Heavily Modified and adapted by: Carlos Lordelo
% Last modified: Jan/2017


if ~isempty(spectrumMagDb_frame)
    freqBin = find(ismember(DFTFrequencies, peakFreq_frame));
    newSpecMagDbPeak_frame = (spectrumMagDb_frame(freqBin)) - 1/8 *((spectrumMagDb_frame(freqBin+1)-spectrumMagDb_frame(freqBin-1)).^2) ./ (spectrumMagDb_frame(freqBin-1)+spectrumMagDb_frame(freqBin+1)-2*spectrumMagDb_frame(freqBin));
    newPeakFreq_frame = freqBin - (1/2)*(spectrumMagDb_frame(freqBin+1)-spectrumMagDb_frame(freqBin-1)) ./ (spectrumMagDb_frame(freqBin-1)+spectrumMagDb_frame(freqBin+1)-2*spectrumMagDb_frame(freqBin));
    newPeakFreq_frame = (newPeakFreq_frame-1).*DFTFrequencies(2);   % Transforming from bin to Hz DFTFrequencies(1) = 0 and DFTFrequencies(2) = Fs/nfft
else
    newPeakFreq_frame = [];
    newSpecMagDbPeak_frame = [];
end
end