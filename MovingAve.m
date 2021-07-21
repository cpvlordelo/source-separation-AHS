function c = MovingAve(y, movL)
% Calculation of moving average of an input vector.
% Dowloaded from web and modified by Zhiyao Duan

% change y to a column vector
if size(y, 1) == 1
    y = y';
end

% force movL to be odd
movL = movL-1+mod(movL,2); 
n = length(y);

% this is the core
c = filter(ones(movL,1)/movL,1,y);

% calculate the ends
cbegin = cumsum(y(1:movL-2));
cbegin = cbegin(1:2:end)./(1:2:(movL-2))';
cend = cumsum(y(n:-1:n-movL+3));
cend = cend(end:-2:1)./(movL-2:-2:1)';
c = [cbegin;c(movL:end);cend];