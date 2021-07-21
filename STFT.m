function [ Yw, t, f ] = STFT( y, wFn, L, S, N, varargin)
%STFT Short-Time Fourier Transform
%
%   Input
%       y - discret signal with length n
%       wFn - window function
%       L - [samples] - size of each frame
%       S - [samples] - hop
%       N - size of FFT
%   Return
%       Yw - STFT{y}

if(S < 1)
    error('Hop in samples must be higher than 0');
end

if rem(N, 2)
    error('Number of FFT points must be multiple of 2');
end

n = length(y);
fs = 0;

%Frames calculation
pre = 0; %Zero padding at begining of signal, better resolution at initial samples.
pos = ceil(1/S*(n-L)); %Zero padding at ending of signal, to cover all n samples.
% if nargout > 1
%     pos = pos + floor((L-1)/S);
%     pre = floor((L-1)/S);
% end

if nargin > 5
 fs = varargin{1};
end
    
M = pre + pos; %Total number of frames

if(isnan(M))
   error('Error getting the number of frames (m)');
end

y = [zeros(pre*S, 1); y; zeros((pos-1)*S + L - n, 1)]; %zero-padding

yw = zeros(L, M);

Yw = zeros(N/2+1, M);

%Calling window function size L
%w = wFn(L, S);
 
%Windowing and FFT
for m=1:M,
    yw(:, m) = y((1:L) + (m-1)*S, 1).* wFn;

    Y = fft(yw(:, m), N); %STFT       
    
    Yw(:, m) = Y(1:N/2+1);
end

% if nargout > 1
%     pad = [pre*S (pos-1)*S + L - n];
% end

if nargout > 2
    % % calculate the time and frequency vectors
    t = (L/2:S:n-L/2-1)/fs;
    f = fs/2*linspace(0,1,N/2+1);
    f = f(:);
end

end
