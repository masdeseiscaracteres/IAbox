clc;
clear;
close all;

K = 4;              % Num users
D=diag(2*ones(K,1));          % Demands matrix
A=ones(K);          % Fully connected system (Xnetwork)
nT=5*ones(K,1);     % Tx antennas
nR=5*ones(K,1);     % Rx antennas

SNR_dB=0:5:30;

%% Find a random IA solution
H0=GenerateChannel(nT,nR,D);
[U, V]=RandomBeamforming(H0,D); %Call the RandomBeamforming algorithm to obtain
H=InverseIA(U,V);
for us=1:K
    H{us,us}=H0{us,us};
end

%% Find ZF decoders
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

%% Find MMSE decoders
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

%% Plot
plot(SNR_dB,SR_ZF,'b-')
hold on;
plot(SNR_dB,SR_MMSE,'r--')
xlabel('SNR [dB]');
ylabel('Sum-rate [bit/s/Hz]');
grid on;
title('ZF vs. MMSE linear decoding')