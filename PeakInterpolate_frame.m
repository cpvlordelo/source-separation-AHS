function [newPeakFreq_frame, newspecMagDb_frame] = PeakInterpolate_frame(specMagDb_frame, peakFreq_frame, DFTFrequencies)
% Quadratic interpolation
%
% Input
%   - specMagDb               : magnitude spectrogram
%   - freqBin                 : frequency bins of detected peaks in a frame
%   - DFTfrequencies          : vector with frequencies (in Hz) corresponding to the
%                               DFT frequency bins
% Output
%   - newPeakFreq_cell        : cell with real frequency (in Hz) of the new detected peaks (each frame is a column)
%   - newspecMagDbPeak_cell   : amplitudes in dB of the new interpolated peaks
% 
% Author: Zhiyao Duan
% Created: 2007
%
% Heavily modified by Carlos Lordelo in 2016


if ~isempty(specMagDb)
    %freqBin = find(ismember(DFTFrequencies, peakFreq));
    newspecMagDb = specMagDb(freqBin) - 1/8 *((specMagDb(freqBin+1)-specMagDb(freqBin-1)).^2) ./ (specMagDb(freqBin-1)+specMagDb(freqBin+1)-2*specMagDb(freqBin));
    newPeakFreq_frame = freqBin - 1/2*(specMagDb(freqBin+1)-specMagDb(freqBin-1)) ./ (specMagDb(freqBin-1)+specMagDb(freqBin+1)-2*specMagDb(freqBin));
else
    newPeakFreq = [];
    newspecMagDb = [];
end
end
