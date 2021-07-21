function [ y ] = ISTFT( Yw, wFn, L, S, zeroPhase, varargin )
%STFT Short-Time Fourier Transform
%
%   Input
%   Yw - STFT{yw}
%   wFn - window function
%   L - [samples] - size of each frame
%   S - [samples] - hop
%   zeroPhase - 
%
%   Return
%   y - ISTFT{Yw}

if(S < 1)
    %error
end

pad = [0 0];

if nargin > 5
    pad = varargin{1};
end

%frames and fft size
[N, M] = size(Yw);

NFFT = (N - 1)*2;

%Total length
n = (M-1)*S + L;

y = zeros(n, 1); %zero-padding x

yw = zeros(NFFT,M);

%Calling window function size L
%w = wFn(L, S, zeroPhase);
    
%Windowing and FFT
    for m=1:M
        
        yw(:, m) = real(ifft([Yw(:, m); conj(Yw(end-1:-1:2, m))], NFFT));
        
        y((1:L) + (m-1)*S) = y((1:L) + (m-1)*S, 1) + yw(1:L, m).*wFn; 
        
    end
    
    w_square_sum = sum(wFn.^2);

    y = y .*S/w_square_sum;

    %Cutting the edges for perfect reconstruction
    y(1:pad(1)) = [];
    y((end-pad(2)+1):end) = [];
end