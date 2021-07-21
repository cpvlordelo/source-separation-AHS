function plotPianoroll(t, F0_hz, parameters )
% Plot the pianoroll of each of the sources 
%
% Input:
%   - t          	     : vector with the time (in seconds) with the center position of each frame of DFT
%   - F0_hz          	 : a matrix where each column has the detected peak's amplitudes of a frame (in dB)
%   - parameters            
%         - Fs           : lowest possible frequency of a F0 (midi number)
%         - frameLen     : highest possible frequency of a F0 (midi number)

Fs = parameters.Fs;
frameLen_bin = parameters.frameLen;

numberSources = size(F0_hz,1);
numberFrames = length(t);

if numberFrames ~= size(F0_hz,2),
    error('The time vector and the F0 vector must have the same length');
end

FrameLen_sec = frameLen_bin./Fs;
t_onset = t(:) - FrameLen_sec/2;               % Do not forget to transform row vectors into column vectors to create nmat
t_duration = FrameLen_sec.*(ones(numberFrames,1));

for i = 1:numberSources,
    idx = find(F0_hz(i,:) > 0);
    pitch_midi = round(hz2midi(F0_hz(i,idx)))';  % Do not forget to transform row vectors into column vectors to create nmat
    numberFrames = length(pitch_midi);
    
%   NMAT FORMAT:  [     ONSET(BEATS)           DURATION (BEATS)          MIDI CHANNEL        MIDI PITCH          VELOCITY             ONSET(SEC)        DURATION(SEC)    ] 
    NMAT =        [zeros(numberFrames, 1),  zeros(numberFrames, 1),  ones(numberFrames, 1),  pitch_midi,  63*ones(numberFrames,1),   t_onset(idx),     t_duration(idx)   ];
    
    figure;
    pianoroll(NMAT,'num', 'sec');
end

