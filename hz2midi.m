function m = hz2midi(hz)
% Hertz to MIDI note number conversion
% 
% Converts frequency values given in Hertz to MIDI note numbers
% Notes are numbered in semitones with middle C being 60. 
% Midi note 69 (A4) has a frequency of 440 hertz
%
% Input arguments: 
%	hz = Frequency values in hertz
%
% Output: 
%	m = MIDI note numbers
%

if isempty(hz), return; end
m=(69+12 * log(abs(hz)/440)/log(2));
