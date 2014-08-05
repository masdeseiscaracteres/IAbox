clc;
clear;
close all;

SNR_dB=0:5:30;
target_SNR_dB=30;

str='(6x6,1)^11';

[nT nR d K]=system_str2vec(str);
D=diag(d);

H=GenerateChannel(nT,nR);

%% Run Incremental-SNR algorithm
[U, V] = IncrementalSNR(H,target_SNR_dB);

%% Run Gauss-Newton algorithm
[U_nwt,V_nwt] = GaussNewtonMinLeakage(H,D);

%% Compute rate
% with ZF decoders
U_ZF=ComputeDecoders(H,V,D,'ZF');
IL_ZF=sum(sum(InterferenceLeakage(H,U_ZF,V,D)));

SR_ZF=zeros(length(SNR_dB),1);
for iSNR=1:length(SNR_dB);
    noise_var=10^(-SNR_dB(iSNR)/10); %Noise variance at each antenna
    
    %Generate noise covariance matrices
    W=cell(K,1);
    for rx=1:K
        W{rx}=noise_var*eye(nR(rx));
    end
    SR_ZF(iSNR)=sum(Rate(H,U_ZF,V,W));
end

% with MMSE decoders
IL_MMSE=zeros(length(SNR_dB),1);
SR_MMSE=zeros(length(SNR_dB),1);
for iSNR=1:length(SNR_dB);
    noise_var=10^(-SNR_dB(iSNR)/10); %Noise variance at each antenna
    
    %Generate noise covariance matrices
    W=cell(K,1);
    for rx=1:K
        W{rx}=noise_var*eye(nR(rx));
    end
    options.W=W;
    U_MMSE=ComputeDecoders(H,V,D,'MMSE',options);
    IL_MMSE(iSNR)=sum(sum(InterferenceLeakage(H,U_MMSE,V,D)));
    SR_MMSE(iSNR)=sum(Rate(H,U_MMSE,V,W));
end

%% Compute rate Gauss-Newton method
IL_nwt=zeros(length(SNR_dB),1);
SR_nwt=zeros(length(SNR_dB),1);
for iSNR=1:length(SNR_dB);
    noise_var=10^(-SNR_dB(iSNR)/10); %Noise variance at each antenna
    
    %Generate noise covariance matrices
    W=cell(K,1);
    for rx=1:K
        W{rx}=noise_var*eye(nR(rx));
    end
    options.W=W;
    IL_nwt(iSNR)=sum(sum(InterferenceLeakage(H,U_nwt,V_nwt,D)));
    SR_nwt(iSNR)=sum(Rate(H,U_nwt,V_nwt,W));
end
%% Plot
plot(SNR_dB,SR_ZF,'b-')
hold on;
plot(SNR_dB,SR_MMSE,'r--')
plot(SNR_dB,SR_nwt,'k-.')
xlabel('SNR [dB]');
ylabel('Sum-rate [bit/s/Hz]');
legend('Incremental SNR (ZF decoders)','Incremental SNR (MMSE decoders)','Gauss-Newton (ZF decoders','Location','Northwest')
grid on;
title('Incremental SNR algorithm vs. Gauss-Newton')

DoF=3*diff(SR_MMSE')./diff(SNR_dB);
fprintf('Estimated DoF Incremental SNR: %.1f\n',DoF(end))

DoF=3*diff(SR_nwt')./diff(SNR_dB);
fprintf('Estimated DoF Gauss-Newton: %.1f\n',DoF(end))