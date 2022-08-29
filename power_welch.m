function spec = power_welch(signal, win_type, win_length, hop_length, K)
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
sig_length = length(signal);
%# reshape the window and signal into a uniform format
window = reshape(window, win_length, 1);
signal = reshape(signal, sig_length, 1);
%# if signal can not be devided by hop length, zero padding
pad_length = mod(sig_length-win_length, hop_length);
if pad_length~=0
    signal = [signal; zeros(hop_length - pad_length, 1)];
    sig_length = length(signal);
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
for i = 1:nums
    x = window.*signal((i-1)*hop_length+1:(i-1)*hop_length+win_length);
    for k = 1:K
        fk = F(:,k);
        estimate_Mat(k, i) = (abs(fk'*x))^2;
    end   
end
%# calculate mean results over all segment estimates
spec_all = mean(estimate_Mat, 2)/win_length;
%# for real signal, one-side spectrum is returned
if isreal(signal)
    spec=spec_all(1:K/2+1);
else
    spec=spec_all;
end
%# recover the scale for different windows
if win_type == "hamming"
    spec=spec*1.586^2; 
elseif win_type == "hanning"
    spec=spec*1.633^2; 
end
%----------------------------------EOF-------------------------------------
