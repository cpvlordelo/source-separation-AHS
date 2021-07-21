function [ peakAmpdB_frameOK, peakFreq_frameOK] = PeakValidation(peakAmpdB_frame, peakFreq_frame, specMagdB_frame, f)
% Validate the peaks found using frame energy
%
% Input
%   - peakAmpdB_frame         : a column vector with the amplitudes (in dB) of the peaks in the frame
%   - peakFreq_frame          : a column vector with the frequencies (in dB) of the peaks in the frame
%   - specmagDb_frame         : frame's magnitude spectrogram (in dB) put in a column vector
%   - f                       : corresponding frequencies (in Hz) of each value in the spectrogram
%
% Output.
%   - peakAmpdB_frameOK       : a column vector with the validated amplitudes (in dB) of the peaks in the frame
%   - peakFreq_frameOK        : a column vector with the validated frequencies (in dB) of the peaks in the frame
%
% Author: Carlos Lordelo
% Last modified: Fev/2017

if length(peakAmpdB_frame) ~= length(peakFreq_frame),
    error('The length of the peaks frequency vector should be the same as the length of the peaks amplitude vector');
else if length(specMagdB_frame) ~= length(f),
    error('The length of the magnitude spectrum should be the same as the length of frequency vector');
    end
end

numberPeaks = numel(peakAmpdB_frame);
peakAmpdB_frameOK = peakAmpdB_frame;
peakFreq_frameOK = peakFreq_frame;
if numberPeaks == 0,
    return;
end
if peakAmpdB_frame(1) <= 0,
    peakAmpdB_frameOK = [];
    peakFreq_frameOK = [];
end
end

