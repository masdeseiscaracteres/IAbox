clc;
clear;
close all;

% K = 6;              % Num users
% T = 1;
% D=diag(5*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=8*ones(K,1);     % Tx antennas
% nR=26*ones(K,1);     % Rx antennas

% K = 4;              % Num users
% T = 1;
% D=diag(6*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=15*ones(K,1);     % Tx antennas
% nR=15*ones(K,1);     % Rx antennas

% K = 5;              % Num users
% T = 1;
% D=diag(10*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=30*ones(K,1);     % Tx antennas
% nR=30*ones(K,1);     % Rx antennas

% K = 3;              % Num users
% T = 1;
% D=diag(2*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=2*ones(K,1);     % Tx antennas
% nR=6*ones(K,1);     % Rx antennas

% K = 3;              % Num users
% T = 1;
% D=diag(1*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=2*ones(K,1);     % Tx antennas
% nR=2*ones(K,1);     % Rx antennas

% K = 4;              % Num users
% T = 1;
% D=diag(8*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=11*ones(K,1);     % Tx antennas
% nR=29*ones(K,1);     % Rx antennas

% K = 3;              % Num users
% T = 1;
% D=diag(3*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=6*ones(K,1);     % Tx antennas
% nR=6*ones(K,1);     % Rx antennas

% K = 3;              % Num users
% T = 1;
% D=diag(4*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=10*ones(K,1);     % Tx antennas
% nR=6*ones(K,1);     % Rx antennas

% K = 19;              % Num users
% T = 1;
% D=diag(1*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=10*ones(K,1);     % Tx antennas
% nR=10*ones(K,1);     % Rx antennas

% K = 3;              % Num users
% T = 1;
% D=diag([3 3 3]);          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=4*ones(K,1);     % Tx antennas
% nR=8*ones(K,1);     % Rx antennas

K = 4;              % Num users
T = 1;
D=diag(2*ones(K,1));          % Demands matrix
A=ones(K);          % Fully connected system (Xnetwork)
nT=5*ones(K,1);     % Tx antennas
nR=5*ones(K,1);     % Rx antennas
 
% K = 5;              % Num users
% T = 1;
% D=diag(1*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=3*ones(K,1);     % Tx antennas
% nR=3*ones(K,1);     % Rx antennas

%% Simulation parameters
SNR_dB=0:20:40;

%%
options.NumExtensions=T;
options.ACS=false;
options.ConstantExtensions=true;
options.A=A;
H = GenerateChannel(nT,nR,options);
[U0,V0] = RandomBeamforming(H,D);

%% MinIL
options.Verbose=0;
options.StartingPoint.V=V0;
options.StartingPoint.U=U0;
options.MaxIter=5000;

[U1, V1]=AlternatingMinLeakage(H,D,options);
 
SR_MinIL=zeros(length(SNR_dB),1);
for iSNR=1:length(SNR_dB);
    noise_var=10^(-SNR_dB(iSNR)/10); %Noise variance at each antenna
    
    %Generate noise covariance matrices
    W=cell(K,1);
    for rx=1:K
        W{rx}=noise_var*eye(nR(rx));
    end    
    SR_MinIL(iSNR)=sum(Rate(H,U1,V1,W));
end
%% MaxSINR
options.Verbose=100;
options.StartingPoint.V=V0;
% options.MaxIter=3000;
options.Orthogonalize=true;

SR_MaxSINR=zeros(length(SNR_dB),1);
for iSNR=1:length(SNR_dB);
    noise_var=10^(-SNR_dB(iSNR)/10); %Noise variance at each antenna
    
    %Generate noise covariance matrices
    W=cell(K,1);
    for rx=1:K
        W{rx}=noise_var*eye(nR(rx));
    end
    options.W=W;
    
    [U2, V2]=AlternatingMaxSINR(H,D,W,options);

    SR_MaxSINR(iSNR)=sum(Rate(H,U2,V2,W));
end

%% MaxSINR_GF
options.Verbose=100;
options.StartingPoint.V=V0;
% options.MaxIter=3000;
options.Orthogonalize=false;

SR_MaxSINR_GF=zeros(length(SNR_dB),1);
for iSNR=1:length(SNR_dB);
    noise_var=10^(-SNR_dB(iSNR)/10); %Noise variance at each antenna
    
    %Generate noise covariance matrices
    W=cell(K,1);
    for rx=1:K
        W{rx}=noise_var*eye(nR(rx));
    end
    options.W=W;
    
    [U3, V3]=AlternatingMaxSINR_GF(H,D,W,options);

    SR_MaxSINR_GF(iSNR)=sum(Rate(H,U3,V3,W));
end

%% Plot
plot(SNR_dB,SR_MinIL,'b-o')
hold on;
plot(SNR_dB,SR_MaxSINR,'r-+')
plot(SNR_dB,SR_MaxSINR_GF,'m-x')
xlabel('SNR [dB]');
ylabel('Sum-rate [bit/s/Hz]');
lh=legend('Min. Leakage','Max. SINR separate filtering','Max. SINR group filtering');
set(lh,'Location','NorthWest');
grid on;
title('Min. Leakage vs. Max. SINR')

% L=1; %Random alignment solution
% SR_RA=AsymptoticSumRateApproxSingleBeam(nT,nR,L,SNR_dB); %random alignment
% plot(SNR_dB,SR_RA,'k:','LineWidth',2);
% L=216; %Max-of-L alignment solution
% SR_MLA=AsymptoticSumRateApproxSingleBeam(nT,nR,L,SNR_dB); %random alignment
% plot(SNR_dB,SR_MLA,'k:','LineWidth',2);
