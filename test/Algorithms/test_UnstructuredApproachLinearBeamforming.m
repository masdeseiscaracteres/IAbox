clc;
clear;
close all;

%% Define scenario
% K = 4;              % Num users
% D=diag(8*ones(K,1));% Demands matrix
% A=ones(K);          % Fully connected system
% nT=11*ones(K,1);     % Tx antennas
% nR=29*ones(K,1);     % Rx antennas

% K = 5;              % Num users
% D=diag(4*ones(K,1));% Demands matrix
% A=ones(K);          % Fully connected system
% nT=19*ones(K,1);     % Tx antennas
% nR=5*ones(K,1);     % Rx antennas

K = 3;              % Num users
D=diag(3*ones(K,1));% Demands matrix
A=ones(K);          % Fully connected system
nT=8*ones(K,1);     % Tx antennas
nR=4*ones(K,1);     % Rx antennas

K = 4;              % Num users
D=diag(3*ones(K,1));% Demands matrix
A=ones(K);          % Fully connected system
nT=4*ones(K,1);     % Tx antennas
nR=11*ones(K,1);     % Rx antennas
% 
% K = 4;              % Num users
% D=diag(1*ones(K,1));% Demands matrix
% A=ones(K);          % Fully connected system
% nT=4*ones(K,1);     % Tx antennas
% nR=1*ones(K,1);     % Rx antennas

%% Generate random channel
H = cell(K,K); %Channel's cell array
for rx = 1:K
    for tx = 1:K
        H{rx,tx}=crandn(nR(rx),nT(tx)); %Random channel matrix
    end
end

%% Solve

[U, V]=UnstructuredApproachLinearBeamforming(H,D);

%% Check residual interference
I = InterferenceLeakage(H,U,V,D)

[r s]=rank_test(D,U,V,H)
s{:}

