clc;
clear;
close all;

%% Does not work in this scenario (0 cost function does not imply 0
%interference leakage)
% K = 4;              % Num users
% D=diag(8*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system
% nT=28*ones(K,1);     % Tx antennas
% nR=12*ones(K,1);     % Rx antennas

%% Gets stagnated
% K = 5;              % Num users
% D=diag(1*ones(K,1));          % Demands matrix
% A=ones(K);          % Fully connected system
% nT=3*ones(K,1);     % Tx antennas
% nR=3*ones(K,1);     % Rx antennas

K = 4;              % Num users
D=diag(1*ones(K,1));          % Demands matrix
A=ones(K);          % Fully connected system
nT=2*ones(K,1);     % Tx antennas
nR=3*ones(K,1);     % Rx antennas
%%
K = 4;              % Num users
D=diag(1*ones(K,1));          % Demands matrix
A=ones(K);          % Fully connected system
nT=3*ones(K,1);     % Tx antennas
nR=2*ones(K,1);     % Rx antennas

%% Generate random channels and initialization point
options.A=A;
H = GenerateChannel(nT,nR,options);
V0=arrayfun(@(a,b)orth(crandn(a,b)),nT',diag(D)','UniformOutput',false);

%%
options.Verbose=1;
options.MaxIterAlternating=500;
options.MaxIterSteepestDescent=50;
options.minUpdate=1e-6;
options.Tol=1e-15;
options.ImprovTol=1e-15;
[U, V, C] = SteepestDescentMinDistance(H,D,options);

I=InterferenceLeakage(H,U,V,D);
IL=sum(sum(I))