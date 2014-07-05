%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     Interference leakage minimization for structured channels       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%% System parameters

K=4; % number of users
D=diag(3*ones(K,1)); % streams per user
nT=1*ones(K,1); % transmit antennas
nR=1*ones(K,1); % receive antennas
L=8; %Number of channel extensions
Improper=false; %If true use assymetric complex signaling
ConstantChannels=false; % If true, constant channel extensions

%%
% K=3; % number of users
% D=diag(1*ones(K,1)); % streams per user
% nT=2*ones(K,1); % transmit antennas
% nR=2*ones(K,1); % receive antennas
% L=1; %Number of channel extensions
% Improper=false; %If true use assymetric complex signaling
% ConstantChannels=false; % If true, constant channel extensions

%% Generate channels
H=generatechannel(nT,nR,ones(K),struct('NumExtensions',L,'ACS',Improper,'ConstantExtensions',ConstantChannels));

%%
options.MaxIter=5000; % number of iterations
options.Tol=1e-8;   % tolerance level for interference leakage
options.L=L;        % number of symbol extensions
options.Improper=Improper; % If true, perform improper signaling with d real streams
options.epsilon=1e-3; % Minimum eigenvalue constraint
options.Verbose=1;

[U,V]=AlternatingMinLeakageStructuredChannel(H,D,options)

%% Sum-rate computation
SNR=0:60;
SR=zeros(size(SNR));
for i=1:length(SNR)
    sigma2=10^(-SNR(i)/10);
    W=cell(1,K);
    for rx=1:K
        W{rx}=sigma2*eye(L*nR(rx));
    end
    SR(i)=sum(Rate(H,U,V,W));
end

%% Plot

if Improper
    figure, plot(SNR,SR/(2*L),'m','linewidth',2), grid on;
    xlabel('SNR [dB]'), ylabel('Sum-Rate [bit/s/Hz]');
    D=diff(SR/(2*L))*3./diff(SNR); % DoF at different SNR values
    figure, plot(SNR,[0 D],'linewidth',2), grid on;
    xlabel('SNR [dB]'), ylabel('Sum-DoF');
else
    figure, plot(SNR,SR/L,'m','linewidth',2), grid on;
    xlabel('SNR [dB]'), ylabel('Sum-Rate [bit/s/Hz]');
    D=diff(SR/L)*3./diff(SNR); % DoF at different SNR values
    figure, plot(SNR,[0 D],'linewidth',2), grid on;
    xlabel('SNR [dB]'), ylabel('Sum-DoF');
end
