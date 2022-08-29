function MSC= coherence_MVDR(signal, param)
%% ------------------------------------------------------------------------
% This project computes the coherence function between two signals with a
% MVDR spectrum estimator
% Usage:
%   MSC = coherence_MVDR(signal, param)
% Ouput:
%   MSC: magnitude squared coherence
% Inputs:
%   signal.Rx: covariance matrix Rx = E(xx')
%   signal.Ry: covariance matrix Ry = E(yy')
%   signal.Rxy: covariance matrix Rxy = E(xy')
%   param.L: length of the first subfilter
%   param.K: fft points
%   param.delta1,delta2: constant for matrix inverse
% Author :
%   Xianrui Wang, Center of Intelligent Acoustics and Immersive
%   Communications(CIAIC)
% Contact:
%   wangxianrui@mail.nwpu.edu.cn
% Reference:
%   ESTIMATION OF THE COHERENCE FUNCTION WITH THE MVDR APPROACH, ICASSP,
%   Jacob Benesty, Jingdong Chen, Arden Huang, 2006.
%--------------------------------------------------------------------------
if nargin<2
    error('Please pass signal and paramter structures');
end
Rx = signal.Rx;           % covariance matrix Rx = E(xx')
Ry = signal.Ry;           % covariance matrix Ry = E(yy')
Rxy = signal.Rxy;         % covariance matrix Rxy = E(xy')
L = param.L;              % length of filter
K = param.K;              % length of frequency bins            
MSC = zeros(1,K);         % magnitude squared coherence function
%--------------------------------------------------------------------------
%% construct fourier matrix
% equation above eq.3, under eq.2
F = zeros(L,K);
l = (0:L-1)';
f = exp(2*pi*l*1j/K);
for k = 0:K-1
    F(:,k+1) = f.^k;
end
F = F/sqrt(L);
%--------------------------------------------------------------------------
%% estimate MSC with MVDR approach
I_Mat = eye(L);             % identity matrix associate with inversion
diag_load_mode = param.diag_load_mode;   % diagonal loading 
if diag_load_mode == "small"             % tr(R)*1e-6
   diag_fac_Rx = trace(Rx)/L*1e-6;
   diag_fac_Ry = trace(Ry)/L*1e-6;
elseif diag_load_mode == "signal_dependent"
   diag_vec_Rx = diag(Rx);                   
   diag_fac_Rx = sqrt(var(diag_vec_Rx));  % std(diag_vec)
   diag_vec_Ry = diag(Ry);                   
   diag_fac_Ry = sqrt(var(diag_vec_Ry));  % std(diag_vec)
else
    error("unsupport diagonal loading way");
end
invRx = (Rx + diag_fac_Rx*I_Mat) \ I_Mat;
invRy = (Ry + diag_fac_Ry*I_Mat) \ I_Mat;
% This part can be replaced by below annotated part for efficiency.
for k = 1:K
    fk = F(:,k);
    %# calculate spectrum of input signals
    % eq.7
    %specx(k+1) = 1/real(fk'*invRx*fk);
    %specy(k+1) = 1/real(fk'*invRy*fk);
    gx = invRx*fk/(fk'*invRx*fk);
    gy = invRy*fk/(fk'*invRy*fk);
    MSC(k) = abs((gx'*Rxy*gy).^2/((gx'*Rx*gx)*(gy'*Ry*gy))); 
end
if isreal(Rxy)
    MSC=MSC(1:K/2+1);
end
%--------------------------------------------------------------------------
%%
%This part is more compuation efficient than line 50-60
% numerator of eq. 20
% numer = ((diag(F'*invRxRxyinvRy*F)));
% denom1 = real(diag(F'*invRx*F));
% denom2 = real(diag(F'*invRy*F));
% % eq. 20
% MSC = real(numer.*conj(numer));%./(denom1.*denom2+denomeps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EOF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
