%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This test script provides usage examples for the Gauss-Newton algorithm %
% in the file GaussNewtonMinLeakage.m                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

%% Some example scenarios
%% (5x5,2)^4
K = 4;               % Number of users
nT=5*ones(K,1);      % Tx antennas
nR=5*ones(K,1);      % Rx antennas
D=diag(2*ones(K,1)); % Demands matrix

K = 3;               % Number of users
nT=2*ones(K,1);      % Tx antennas
nR=2*ones(K,1);      % Rx antennas
D=diag(1*ones(K,1)); % Demands matrix
%% (3x3,1)^5
% K = 5;
% nT=3*ones(K,1);
% nR=3*ones(K,1);
% D=diag(1*ones(K,1));
%% (10x10,1)^19
% K = 19;
% nT=10*ones(K,1);
% nR=10*ones(K,1);
% D=diag(1*ones(K,1));
%% (12x12,4)^5
% K = 5;
% nT=12*ones(K,1);
% nR=12*ones(K,1);
% D=diag(4*ones(K,1));
%% (29x11,8)^4
% K = 4;
% nT=29*ones(K,1);
% nR=11*ones(K,1);
% D=diag(8*ones(K,1));
%% (2x3,1)^2(3x2,1)^2
% K=4;
% nT=[3 3 2 2];
% nR=[2 2 3 3];
% D=diag([1 1 1 1]);

%%
K=3;
nT=[9 9 9];
nR=nT;
D=3*ones(K,1);
%% Prepare input arguments
% Generate a random interference channel
H = cell(K,K); %Channel's cell array
for rx = 1:K
    for tx = 1:K
        H{rx,tx}=crandn(nR(rx),nT(tx)); %Random channel matrix
    end
end
% Define some options
options.NwtTol=1e-15;
options.Verbose=0;

%% Now run the algorithm
tic;
[U,V,IL,success] = GaussNewtonMinLeakage(H,D,options);
t=toc;
fprintf('Success: %d, Elapsed time: %.2f, Interference leakage: %.2e\n',success,t,IL(end));

%% Plot convergence curve
semilogy(0:length(IL)-1,IL,'b');
grid on;
xlabel('Iteration number');
ylabel('Interference leakage');
title('Interference leakage evolution with the number of Gauss-Newton iterations');
