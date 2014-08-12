clc;
clear;
close all;

%% Some example scenarios
% K = 4;              % Num users
% D=diag(8*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system
% nT=28*ones(K,1);     % Tx antennas
% nR=12*ones(K,1);     % Rx antennas

K = 5;              % Num users
D=diag(1*ones(K,1));          % Demands matrix
A=ones(K);          % Fully connected system
nT=3*ones(K,1);     % Tx antennas
nR=3*ones(K,1);     % Rx antennas

% K = 4;              % Num users
% D=diag(1*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system
% nT=2*ones(K,1);     % Tx antennas
% nR=3*ones(K,1);     % Rx antennas

% K = 4;              % Num users
% D=diag(1*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system
% nT=3*ones(K,1);     % Tx antennas
% nR=2*ones(K,1);     % Rx antennas

% K = 4;              % Num users
% D=diag(2*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system
% nT=5*ones(K,1);     % Tx antennas
% nR=5*ones(K,1);     % Rx antennas

%% Generate random channels and maximum power initialization point
options.NumExtensions=1; %No channel extensions
options.ACS=false;
options.ConstantExtensions=true;
options.A=A;
H = GenerateChannel(nT,nR,options);
[U,V]=RandomBeamforming(H,D);

% V=cell(K,1);
% for tx=1:K
%   [F,S,G]=svd(H{tx,tx});
%   V{tx}=G(:,1:D(tx,tx));
% end
%% Run algorithm

options.MaxIter=5000;
options.Tstart=1;
options.tol=1e-10;
options.eta=0.995;
options.Verbose=0;
options.StartingPoint.U=U;
options.StartingPoint.V=V;
[U,V,Power]=AlternatingMaxPower(H,D,options);

[U2,V2,Power2]=AlternatingMinLeakage(H,D,options);

%Equivalently
% options.eta=0;
% [U2,V2,Power3]=AlternatingMaxPower(H,D,options);

%% Calculate sum-rate
SNR_dB=0:5:40;
SR=zeros(length(SNR_dB),1);
SR2=zeros(length(SNR_dB),1);
for iSNR=1:length(SNR_dB);
    noise_var=10^(-SNR_dB(iSNR)/10); %Noise variance at each antenna
    
    %Generate noise covariance matrices
    W=cell(K,1);
    for rx=1:K
        W{rx}=noise_var*eye(nR(rx));
    end
    options.W=W;
    
    SR(iSNR)=sum(Rate(H,U,V,W));
    SR2(iSNR)=sum(Rate(H,U2,V2,W));
end

%% Plot
plot(SNR_dB,SR,'b-o')
hold on;
plot(SNR_dB,SR2,'r-+')
xlabel('SNR [dB]');
ylabel('Sum-rate [bit/s/Hz]');
lh=legend('Max. power','Min. Leakage');
set(lh,'Location','NorthWest');
grid on;
title('Sum-rate performance of the Max. Power solution')
