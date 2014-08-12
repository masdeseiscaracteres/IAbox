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

K = 3;              % Num users
T = 1;
D=diag(3*ones(K,1));          % Demands matrix
A=ones(K);          % Fully connected system (Xnetwork)
nT=6*ones(K,1);     % Tx antennas
nR=6*ones(K,1);     % Rx antennas

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

% K = 3;              % Num users
% T = 5;              % Time extensions
% D=diag(4*ones(K,1));  % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=1*ones(K,1);     % Tx antennas
% nR=1*ones(K,1);     % Rx antennas
%
% K = 4;              % Num users
% T = 3;              % Time extensions
% D=diag([2 2 2 2]);  % Demands matrix
% A=ones(K);          % Fully connected system (Xnetwork)
% nT=1*ones(K,1);     % Tx antennas
% nR=1*ones(K,1);     % Rx antennas

%% Simulation parameters
SNR_dB=-20:5:60;

%%
options.NumExtensions=T;
options.ACS=false;
options.ConstantExtensions=true;
options.A=A;
H = GenerateChannel(nT,nR,options);
[U0,V0] = RandomBeamforming(H,D);

options.Verbose=0;
options.StartingPoint.U=U0;
options.StartingPoint.V=V0;
options.MaxIter=5000;
[U1, V1, IL1]=AlternatingMinLeakage(H,D,options);

options.Verbose=0;
options.StartingPoint.U=U0;
options.StartingPoint.V=V0;
options.w=0.01;
options.EndingUpdate=1e-4;
options.MaxIter=5000;
[U2, V2, IL2]=HybridSignalInterferenceAligment(H,D,options);
%% Compute sum-rate performance
SR1=zeros(length(SNR_dB),1);
SR2=zeros(length(SNR_dB),1);
for iSNR=1:length(SNR_dB);
    noise_var=10^(-SNR_dB(iSNR)/10); %Noise variance at each antenna
    
    %Generate noise covariance matrices
    W=cell(K,1);
    for rx=1:K
        W{rx}=noise_var*eye(nR(rx));
    end
    options.W=W;
    SR1(iSNR)=sum(Rate(H,U1,V1,W));
    SR2(iSNR)=sum(Rate(H,U2,V2,W));
end
%% Plot
plot(SNR_dB,SR1,'b-')
hold on;
plot(SNR_dB,SR2,'r--')
xlabel('SNR [dB]');
ylabel('Sum-rate [bit/s/Hz]');
lh=legend('Min. Leakage','Joint Signal+Interference optim.');
set(lh,'Location','NorthWest');
grid on;
title('MinLeakage vs. Joint Signal+Interf. optim.')