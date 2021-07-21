function [SDR, SIR, SAR] = Evaluate(OriginalSourceFile, SeparatedSourceFile, varargin)

addpath ..\..\Codes\Matlab\BSS_EvalToolbox;
s = audioread(OriginalSourceFile);
x = audioread(SeparatedSourceFile);

s = s(:)'; % Making sure s is a row vector
x = x(:)'; % Making sure x is a row vector

if nargin > 2, % We should ignore some parts at the beginning and at the ending of the signals
    s(1:varargin{1}) = [];
    x(1:varargin{1}) = [];

    s(end:-1:end - varargin{1} +1) = [];
    x(end:-1:end - varargin{1} +1) = [];
end

%% Putting the length of the signals the same
if length(x) > length(s),
    x = x(1:length(s));
else
    s = s(1:length(x));
end

%% Evaluating   
%S = [organ'; piccolo'];

[s_target, e_interf, e_artif] = bss_decomp_gain(x, 1, s);
[SDR, SIR, SAR] = bss_crit(s_target, e_interf, e_artif);

end