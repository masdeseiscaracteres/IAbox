%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     Interference leakage minimization for structured channels       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% System parameters
% Clarifying note:
% If Improper=true, the demands matrix is considered to contain the 
% number of real streams instead of complex streams

%% Examples taken from:
% C. Lameiro, Ó. González, and I. Santamaría, "An Interference Alignment 
% Algorithm for Structured Channels," in IEEE 14th Workshop on Signal 
% Processing Advances in Wireless Communications (SPAWC), 2013, pp. 295?299.

% [(1x1,3)^4,8], proper signaling, time-varying channel extensions, 1.5 DoF
K=4; % Number of users
D=diag(3*ones(K,1)); % Streams per user
nT=1*ones(K,1); % Transmit antennas
nR=1*ones(K,1); % Receive antennas
L=8; % Number of channel extensions
Improper=false; % If true use assymetric complex signaling
ConstantChannels=false; % If true, constant channel extensions

% [(1x1,1)^4,3], improper signaling, constant channels, 1.33 DoF
% K=4;
% D=diag(2*ones(K,1));
% nT=1*ones(K,1);
% nR=1*ones(K,1);
% L=3;
% Improper=true;
% ConstantChannels=true;

% [(2x1,2)^3(2x1,3)^3,6], proper signaling, time-varying channels, 2.5 DoF
% K=6;
% D=diag([2 2 2 3 3 3]);
% nT=2*ones(K,1);
% nR=1*ones(K,1);
% L=6;
% Improper=false;
% ConstantChannels=false;

%% More examples

% K=3;
% D=diag(1*ones(K,1));
% nT=2*ones(K,1);
% nR=2*ones(K,1);
% L=1;
% Improper=false;
% ConstantChannels=false;

%% Generate channels
H=GenerateChannel(nT,nR,ones(K),struct('NumExtensions',L,'ACS',Improper,'ConstantExtensions',ConstantChannels));

%% Run the algorithm
options.MaxIter=500; % number of iterations
options.Tol=1e-8;   % tolerance level for interference leakage
options.L=L;        % number of symbol extensions
options.Improper=Improper; % If true, perform improper signaling with d real streams
options.epsilon=1e-3; % Minimum eigenvalue constraint
options.Verbose=100;

F=1+options.Improper;
fprintf('Requested sum-DoF per channel use: %.2f\n',sum(D(:))/L/F);
[U,V]=AlternatingMinLeakageStructuredChannel(H,D,options)

%% Sum-rate computation
SNR=0:60;
SR=zeros(size(SNR));
for i=1:length(SNR)
    sigma2=10^(-SNR(i)/10);
    W=cell(1,K);
    for rx=1:K
        W{rx}=sigma2*eye((1+options.Improper)*L*nR(rx));
    end
    SR(i)=sum(Rate(H,U,V,W));
end

%% Plot

figure, plot(SNR,SR/(F*L),'m','linewidth',2), grid on;
xlabel('SNR [dB]'), ylabel('Sum-Rate [bit/s/Hz]');
D=diff(SR/(F*L))*3./diff(SNR); % DoF at different SNR values
figure, plot(SNR(2:end),D,'linewidth',2), grid on;
xlabel('SNR [dB]'), ylabel('Complex Sum-DoF');

