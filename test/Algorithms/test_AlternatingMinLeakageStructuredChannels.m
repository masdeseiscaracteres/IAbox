%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     Interference leakage minimization for structured channels       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Examples:
% Examples taken from:
% C. Lameiro, Ó. González, and I. Santamaría, "An Interference Alignment 
% Algorithm for Structured Channels," in IEEE 14th Workshop on Signal 
% Processing Advances in Wireless Communications (SPAWC), 2013, pp. 295?299.

% [nT, nR, D, options]=GetTestScenario('[(1x1,3)^4,8]');
% [nT, nR, D, options]=GetTestScenario('[(1x1,1)^4,3]');
[nT, nR, D, options]=GetTestScenario('[(2x1,2)^3(2x1,3)^3,6]');

% Clarifying note:
% If Improper=true, the demands matrix is considered to contain the 
% number of real streams instead of complex streams

%Define some auxiliary variables
F=1+options.ACS;
L=options.NumExtensions; %Number of channel uses
K=length(nT);

%% Generate channels
H=GenerateChannel(nT,nR,options);

%% Run the algorithm
% Set algorithm parameters
options.MaxIter=5000; % number of iterations
options.Tol=1e-8;   % tolerance level for interference leakage
options.epsilon=1e-3; % Minimum eigenvalue constraint
options.Verbose=100;

fprintf('Requested sum-DoF per channel use: %.2f\n',sum(D(:))/L/F);
[U, V] = AlternatingMinLeakageStructuredChannel(H,D,options);

%% Sum-rate computation
SNR=0:60;
SR=zeros(size(SNR));
for i=1:length(SNR)
    sigma2=10^(-SNR(i)/10);
    W=cell(1,K);
    for rx=1:K
        W{rx}=sigma2*eye(F*L*nR(rx));
    end
    SR(i)=sum(Rate(H,U,V,W));
end

%% Plot

figure, plot(SNR,SR/(F*L),'m','linewidth',2), grid on;
xlabel('SNR [dB]'), ylabel('Sum-Rate [bit/s/Hz]');
D=diff(SR/(F*L))*3./diff(SNR); % DoF at different SNR values
figure, plot(SNR(2:end),D,'linewidth',2), grid on;
xlabel('SNR [dB]'), ylabel('Complex Sum-DoF');

