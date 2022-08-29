function MSC = coherence_welch(x, y, win_type, win_length, hop_length, K)
%% ---------------------------------------------------------------
% This project computes the power spectrum
% Usage:
%   spec = power_welch(x,filter_length, hop_length)
% Ouput:
%   spec: power spectrum
% Inputs:
%   x: input signal
%   win_type: type of window
%   win_length: length of window
%   hop_length: overlap size
%   K: fft points
% Author :
%   Xianrui Wang, Center of Intelligent Acoustics and Immersive
%   Communications(CIAIC)
% Contact:
%   wangxianrui@mail.nwpu.edu.cn
%--------------------------------------------------------------------------
if nargin<5
    error('Please pass all parameters');
end
if win_type == "hamming"
    window = hamming(win_length);
elseif win_type == "rect"
    window = ones(win_length,1);
elseif win_type == "hanning"
    window = hanning(win_length);
else
    error("unsupported window type");
end
sig_length = length(x);
%# reshape the window and signal into a uniform format
window = reshape(window, win_length, 1);
x = reshape(x, sig_length, 1);
%# if signal can not be devided by hop length, zero padding
pad_length = mod(sig_length-win_length, hop_length);
if pad_length~=0
    x = [x; zeros(hop_length - pad_length, 1)];
    y = [y; zeros(hop_length - pad_length, 1)];
    sig_length = length(x);
end
%--------------------------------------------------------------------------
%# partion signal into overlapped segments 
nums = floor((sig_length - win_length)/hop_length) + 1;
estimate_Mat = zeros(K,nums);
%# construct Fourier matrix
F = zeros(win_length,K);           % initial Fourier Matrix
l = (0:win_length-1)';             % length of fft vector
f = exp(2*pi*l*-1j/K);             % fft vector of base frequency
for k = 0:K-1
    % fft vector
    F(:,k+1) = f.^k/sqrt(win_length);      
end
%# estimate power spectrum of every signal segments
Pxx = zeros(K,1); 
Pyy = zeros(K,1); 
Pxy = zeros(K,1); 
for i = 1:nums
    x_vec = window.*x((i-1)*hop_length+1:(i-1)*hop_length+win_length);
    y_vec = window.*y((i-1)*hop_length+1:(i-1)*hop_length+win_length);
    Xx = fft(x_vec, K);
    Yy = fft(y_vec, K);
    Xx2 = abs(Xx).^2;
    Yy2 = abs(Yy).^2;
    Xy2 = Yy.*conj(Xx);
    Pxx = Pxx + Xx2;
    Pyy = Pyy + Yy2;
    Pxy = Pxy + Xy2;
end
MSC = (abs(Pxy)).^2./(Pxx.*Pyy);   

if isreal(x)
    MSC = MSC(1:K/2+1);
end
%----------------------------------EOF-------------------------------------
