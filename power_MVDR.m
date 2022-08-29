function spec = power_MVDR(signal, param)
%% ---------------------------------------------------------------
% This project computes the power spectrum
% Usage:
%   spec = power_MVDR(signal, param)
% Ouput:
%   spec: power spectrum
% Inputs:
%   signal.R: covariance matrix R = E(xx')
%   param.L: length of the first subfilter
%   param.K: fft points
%   param.diag_load_mode: two diagonal loading ways
% Reference:
%   ESTIMATION OF THE COHERENCE FUNCTION WITH THE MVDR APPROACH, ICASSP,
%   Jacob Benesty, Jingdong Chen, Arden Huang, 2006.
% Author :
%   Xianrui Wang, Center of Intelligent Acoustics and Immersive
%   Communications(CIAIC)
% Contact:
%   wangxianrui@mail.nwpu.edu.cn
%--------------------------------------------------------------------------
if nargin<2
    error('Please pass signal and paramter structures');
end
R = signal.R;                            % covariance matrix Rx = E(xx')
L = param.L;                             % length of filter
K = param.K;                             % number of frequency bins 
spec_all = zeros(1,K);                   % two-sides spectrum of signal x
diag_load_mode = param.diag_load_mode;   % diagonal loading 
if diag_load_mode == "small"             % tr(R)*1e-6
   diag_load_fac = trace(R)/L*1e-6;
elseif diag_load_mode == "signal_dependent"
   diag_vec = diag(R);                   
   diag_load_fac = sqrt(var(diag_vec));  % std(diag_vec)
else
    error("unsupport diagonal loading way");
end
%--------------------------------------------------------------------------
%% construct fourier matrix
%# equation above eq.3, under eq.2
F = zeros(L,K);
l = (0:L-1)';
f = exp(2*pi*l*1j/K);
for k = 0:K-1
    F(:,k+1) = f.^k/sqrt(L);      % fft vector
end
%--------------------------------------------------------------------------
%% estimate power spectrum  
I_Mat = eye(L);
invR = (R+diag_load_fac*I_Mat)\I_Mat;
for k = 1:K
    fk = F(:,k);
    %# eq.7, calculate the MVDR filter
    g = (invR*fk)/(fk'*invR*fk);
    %# calculate spectrum
    spec_all(k) = abs(g'*R*g);
end
% recover the scale
spec_all = spec_all/L;
% for real signal, one-side spectrum is returned
if isreal(R)
    spec=spec_all(1:K/2+1);
else
    spec=spec_all;
end
%---------------------------------EOF--------------------------------------