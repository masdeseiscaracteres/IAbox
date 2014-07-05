function SR=AsymptoticSumRateApproxSingleBeam(nT,nR,L,SNR)
% Computes an asymptotic sum-rate approximation for the best out of L
% interference aligment solutions in a single-beam interference channel and
% the SNRs given in vector SNR:
%
% Based on the results in
%
% D. A. Schmidt, W. Utschick, and M. L. Honig, "Large System Performance of
% Interference Alignment in Single-Beam MIMO Networks," in 2010 IEEE Global
% Telecommunications Conference GLOBECOM".

% To Do: generalize to asymmetric scenarios

if ~all(nT==nT(1)) || ~all(nR==nR(1))
   error('This function accepts symmetric interference channels only') 
elseif length(nT)~=length(nR)
   error('The number of transmit/receive users does not agree')  
end

K=length(nT);
nT=nT(1);
nR=nR(1);

D=K; %one stream per user
slope=D/3; %only when SNR in dB, sum-rate in bit/s

gamma=0.5772156649; %Euler-Mascheroni constant

if L==1 %Random aligment solution
    improv=0;
else
    ell=qfuncinv(1./L);
    improv=pi*sqrt((nT+nR-1)/6).*(ell+gamma./ell);
end

r=-(nT+nR-1)*gamma+improv; %Rate offset at SNR=0 dB

SR=r+SNR*slope;