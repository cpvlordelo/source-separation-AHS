function [wavData, fs] = ReadWavFile_Mono(wavFile, TimeLength)
% Read wavfile into a single track, truncate, zero mean, and normalize to maximum value of 1.
% Average the two channels if the file is stereo.
% Input:
%   - wavFile   : file name
%   - TimeLength: time length in ms
% Output:
%   - wavData   : mono data, column vector
%   - fs        : sampling rate
% Author: Zhiyao Duan
% Modified: 5/24/2013

% Modified by Carlos Lordelo in 7/11/2016

[wavData, fs] = audioread(wavFile);
wavData = mean(wavData, 2);         % convert stereo to mono
L = length(wavData);

if nargin == 1,
    wavData = wavData - mean(wavData);  % zero mean
    wavData = wavData./sqrt(mean(wavData.^2));  % normalization RMS
    %wavData = wavData./max(abs(wavData)); % normalization of maximum value of 1
end
if nargin == 2,
    tmpwavLength = floor(TimeLength*fs/1000);
    L = min(L, tmpwavLength);       % truncate wavData if it is longer than input
    wavData = wavData(1:L);
    wavData = wavData - mean(wavData);  % zero mean
    wavData = wavData./sqrt(mean(wavData.^2));  % normalization RMS
    %wavData = wavData./max(abs(wavData)); % normalization of maximum value of 1
end

end



% figure; plot(1:L, wavData);
