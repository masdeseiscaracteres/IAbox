function SR=AverageSumRateApprox(K,d,SNR_dB)
% Approximation from:
% O. El Ayach, A. Lozano, and R. W. Heath, "On the Overhead of 
% Interference Alignment: Training, Feedback, and Cooperation," 
% IEEE Trans. Wireless Commun., vol. 11, no. 11, pp. 4192–4203.
SR=zeros(1,length(SNR_dB));
for iSNR=1:length(SNR_dB)
    SNR=10^(SNR_dB(iSNR)/10);
    eta=1/SNR;
    f=@(t)exp(-eta*t)./t;
    E1=quad(f,1,1e10);
    SR(iSNR)=K*d*log2(exp(1))*exp(1/SNR)*E1;
end