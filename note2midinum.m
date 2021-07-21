function MIDInum = note2midinum(name)
% Convert note name to MIDI number. 
%
% Input:
%   - name      : note name e.g. 'C4' or 'D5#'. You should follow this format.
% Output:
%   - MIDInum   : MIDI number of the note name
%
% Author: Zhiyao Duan
% Created: 2007
% Last modified: 6/26/2012

SemiNote1 = ('CCDDEFFGGAAB')';
SemiNote2 = (' # #  # # # ')';
SemiNote3 = ('CDDEEFGGAABB')';
SemiNote4 = (' b b  b b b ')';

if length(name)<3                                   % e.g. 'C4'
    MIDInum = min(strfind(SemiNote1', name(1)));    % min is the pitch class
elseif name(3) == '#'                               % e.g. 'C4#'
    MIDInum = max(strfind(SemiNote1', name(1)));    % max is the pitch class
elseif name(3) == 'b'                               % e.g. 'D4b'
    MIDInum = min(strfind(SemiNote3', name(1)));    % min is the pitch class
else                                                % treat other cases the same as the first case
    MIDInum = min(strfind(SemiNote1', name(1))); 
end
% 69 = 'A4'
MIDInum = MIDInum-10 + 69 + (str2double(name(2))-4)*12;
